#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "index.h"
#include "util.h"
#include "kseq.h"
#include "zran.h"
#include "sequence.h"

/*
create an index
@param file_path, fasta path and name
@param uppercase, uppercase sequence
@param uppercase
*/
pyfastx_Index* pyfastx_init_index(PyObject *obj, PyObject* file_obj, PyObject* index_obj, int uppercase, int full_name, int memory_index, PyObject* key_func){
	pyfastx_Index* index;

	char *index_file;
	Py_ssize_t index_len;

	index = (pyfastx_Index *)malloc(sizeof(pyfastx_Index));
	index->uppercase = uppercase;

	//key function
	if (key_func) {
		index->key_func = Py_NewRef(key_func);
	} else {
		index->key_func = NULL;
	}

	//full name
	index->full_name = full_name;

	//check input file is gzip or not
	index->gzip_format = is_gzip_format(file_obj);

	//initial kseqs
	index->gzfd = pyfastx_gzip_open(file_obj, "rb");
	index->kseqs = kseq_init(index->gzfd);

	//create index file or memory index
	if (memory_index) {
		index->index_file = (char *)malloc(sizeof(char)*9);
		strcpy(index->index_file, ":memory:");
	} else {
		if (index_obj) {
			index_file = (char *)PyUnicode_AsUTF8AndSize(index_obj, &index_len);
			index->index_file = (char *)malloc(index_len);
			memcpy(index->index_file, index_file, index_len);
			index->index_file[index_len] = '\0';
		} else {
			index_file = (char *)PyUnicode_AsUTF8AndSize(file_obj, &index_len);
			index_len += 5;
			index->index_file = (char *)malloc(index_len);
			strcpy(index->index_file, index_file);
			strcat(index->index_file, ".fxi");
		}
	}

	index->fd = _Py_fopen_obj(file_obj, "rb");

	index->index_db = 0;

	if(index->gzip_format){
		index->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
		//initial zran index
		zran_init(index->gzip_index, index->fd, NULL, 1048576, 32768, 16384, ZRAN_AUTO_BUILD);
	}

	//cache name
	index->cache_chrom = 0;

	//cache start and end position
	index->cache_start = 0;
	index->cache_end = 0;

	//complete seq or not
	index->cache_full = 0;

	//enter iteration loop
	index->iterating = 0;

	//prepared sql
	index->iter_stmt = NULL;
	index->uid_stmt = NULL;
	index->seq_stmt = NULL;

	//cache sequence
	//index->cache_name = {0,0,0};
	//index->cache_seq = {0,0,0};
	kstring_init(index->cache_name);
	kstring_init(index->cache_seq);

	//parent fasta
	index->fasta = obj;

	return index;
}

void pyfastx_rewind_index(pyfastx_Index *self){
	kseq_rewind(self->kseqs);
	gzrewind(self->gzfd);
}

void pyfastx_create_index(pyfastx_Index *self){
	// seqlite3 return value
	int ret;
	
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	
	// 1: normal fasta sequence with the same length in line
	// 0: not normal fasta sequence with different length in line
	int seq_normal = 1;

	//1: \n, 2: \r\n
	int line_end = 1;
	
	//length of previous line
	Py_ssize_t line_len = 0;

	//length of current line
	Py_ssize_t temp_len = 0;

	//number of lines that line_len not equal to temp_len
	Py_ssize_t bad_line = 0;

	//current read position
	Py_ssize_t position = 0;
	
	// start position
	Py_ssize_t start = 0;

	//sequence length
	Py_ssize_t seq_len = 0;

	//real line len
	Py_ssize_t real_len;

	//total sequence count
	Py_ssize_t total_seq = 0;

	//total sequence length
	Py_ssize_t total_len = 0;

	//space position
	//char *space_pos;

	//reading file for kseq
	kstream_t* ks;

	//read for line
	kstring_t line = {0, 0, 0};

	char *header_pos;

	//description line length
	int desc_len = 0;

	//chromosome name
	kstring_t chrom = {0, 0, 0};
	char *temp_chrom;

	const char *sql;

	PYFASTX_SQLITE_CALL(ret = sqlite3_open(self->index_file, &self->index_db));
	
	if (ret != SQLITE_OK) {
		PyErr_Format(PyExc_ConnectionError, "Could not open index file %s", self->index_file);
		return;
	}
	
	//create index database
	sql = " \
		CREATE TABLE seq ( \
			ID INTEGER PRIMARY KEY, --seq identifier\n \
			chrom TEXT, --seq name\n \
			boff INTEGER, --seq offset start\n \
			blen INTEGER, --seq byte length\n \
			slen INTEGER, --seq length\n \
			llen INTEGER, --line length\n \
			elen INTEGER, --end length\n \
			norm INTEGER, --line with the same length or not\n \
			dlen INTEGER --description header line length\n \
		); \
		CREATE TABLE stat ( \
			seqnum INTEGER, --total seq counts \n \
			seqlen INTEGER, --total seq length \n \
			avglen REAL, --average seq length \n \
			medlen REAL, --median seq length \n \
			n50 INTEGER, --N50 seq length \n \
			l50 INTEGER --N50 seq count \n \
		); \
		CREATE TABLE comp ( \
			ID INTEGER PRIMARY KEY, --comp identifier\n \
			seqid INTEGER, --seq id\n \
			abc INTEGER, --seq letter\n \
			num INTEGER -- letter count\n \
		); \
		CREATE TABLE gzindex ( \
			ID INTEGER PRIMARY KEY, \
			content BLOB \
		);";

	PYFASTX_SQLITE_CALL(ret = sqlite3_exec(self->index_db, sql, NULL, NULL, NULL));

	if (ret != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, "could not create index tables");
		return;
	}
	
	sql = "PRAGMA synchronous=OFF; PRAGMA locking_mode=EXCLUSIVE; BEGIN TRANSACTION;";
	PYFASTX_SQLITE_CALL(ret = sqlite3_exec(self->index_db, sql, NULL, NULL, NULL));
	if (ret != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, "could not begin transaction");
		return;
	}

	sql = "INSERT INTO seq VALUES (?,?,?,?,?,?,?,?,?);";
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));
	
	gzrewind(self->gzfd);
	ks = ks_init(self->gzfd);

	//Py_BEGIN_ALLOW_THREADS
	while (ks_getuntil(ks, '\n', &line, 0) >= 0) {
		position += line.l + 1;

		//first char is >
		if (line.s[0] == 62) {
			if (start > 0) {
				//end of sequence and check whether normal fasta
				seq_normal = (bad_line > 1) ? 0 : 1;
				
				PYFASTX_SQLITE_CALL(
					sqlite3_bind_null(stmt, 1);
					sqlite3_bind_text(stmt, 2, chrom.s, chrom.l, SQLITE_STATIC);
					sqlite3_bind_int64(stmt, 3, start);
					sqlite3_bind_int(stmt, 4, position-start-line.l-1);
					sqlite3_bind_int64(stmt, 5, seq_len);
					sqlite3_bind_int64(stmt, 6, line_len);
					sqlite3_bind_int(stmt, 7, line_end);
					sqlite3_bind_int(stmt, 8, seq_normal);
					sqlite3_bind_int(stmt, 9, desc_len);
					sqlite3_step(stmt);
					sqlite3_reset(stmt);
				);

				++total_seq;
				total_len += seq_len;
			}

			//reset
			start = position;
			seq_len = 0;
			temp_len = 0;
			line_len = 0;
			line_end = 1;
			bad_line = 0;
			seq_normal = 1;

			//get line end length \r\n or \n
			if (line.s[line.l-1] == '\r') {
				line_end = 2;
			}

			desc_len = line.l - line_end;

			if (chrom.m < line.l) {
				chrom.m = line.l;
				chrom.s = (char *)realloc(chrom.s, chrom.m);
			}

			//remove > sign at the start position
			header_pos = line.s + 1;

			if (self->key_func == NULL) {
				if (self->full_name) {
					chrom.l = desc_len;
					memcpy(chrom.s, header_pos, chrom.l);
					chrom.s[chrom.l] = '\0';
				} else {
					//space_pos = strchr(header_pos, ' ');
					//find space or tab
					for (chrom.l = 0; chrom.l < desc_len; ++chrom.l) {
						if ((header_pos[chrom.l] == ' ') || (header_pos[chrom.l] == '\t')) {
							break;
						}
					}

					//if (space_pos == NULL) {
					//	chrom.l = desc_len;
					//} else {
					//	chrom.l = space_pos - header_pos;
					//}
					memcpy(chrom.s, header_pos, chrom.l);
					chrom.s[chrom.l] = '\0';
				}
			} else {
				// if in Py_BEGIN_ALLOW_THREADS, first acquire the GIL
				//PyGILState_STATE state = PyGILState_Ensure();
				PyObject *result = PyObject_CallFunction(self->key_func, "s", header_pos);
				//PyGILState_Release(state);

				if (!result) {
					PyErr_Print();
					return;
				}

				temp_chrom = (char *)PyUnicode_AsUTF8AndSize(result, &chrom.l);
				memcpy(chrom.s, temp_chrom, chrom.l);
				chrom.s[chrom.l] = '\0';
				Py_XDECREF(result);
			}

			continue;
		}

		temp_len = line.l + 1;

		if (line_len > 0 && line_len != temp_len) {
			bad_line++;
		}

		//record first line length
		if (line_len == 0) {
			line_len = temp_len;
		}

		//calculate atgc counts
		real_len = line.l - line_end + 1;

		//calculate seq len
		seq_len += real_len;
	}

	//end of sequence and check whether normal fasta
	seq_normal = (bad_line > 1) ? 0 : 1;
	
	PYFASTX_SQLITE_CALL(
		sqlite3_bind_null(stmt, 1);
		sqlite3_bind_text(stmt, 2, chrom.s, chrom.l, SQLITE_STATIC);
		sqlite3_bind_int64(stmt, 3, start);
		sqlite3_bind_int(stmt, 4, position-start);
		sqlite3_bind_int64(stmt, 5, seq_len);
		sqlite3_bind_int64(stmt, 6, line_len);
		sqlite3_bind_int(stmt, 7, line_end);
		sqlite3_bind_int(stmt, 8, seq_normal);
		sqlite3_bind_int(stmt, 9, desc_len);
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	stmt = NULL;

	++total_seq;
	total_len += seq_len;
	
	PYFASTX_SQLITE_CALL(
		sqlite3_exec(self->index_db, "PRAGMA locking_mode=NORMAL;", NULL, NULL, NULL);
		sqlite3_exec(self->index_db, "COMMIT;", NULL, NULL, NULL);
		sqlite3_exec(self->index_db, "CREATE UNIQUE INDEX chromidx ON seq (chrom);", NULL, NULL, NULL);
		sqlite3_prepare_v2(self->index_db, "INSERT INTO stat (seqnum,seqlen) VALUES (?,?);", -1, &stmt, NULL);
		sqlite3_bind_int64(stmt, 1, total_seq);
		sqlite3_bind_int64(stmt, 2, total_len);
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);
	
	//Py_END_ALLOW_THREADS

	ks_destroy(ks);
	free(line.s);
	free(chrom.s);

	//create gzip random access index
	if (self->gzip_format) {
		if (strcmp(self->index_file, ":memory:") == 0) {
			zran_build_index(self->gzip_index, 0, 0);
		} else {
			pyfastx_build_gzip_index(self->gzip_index, self->index_db);
		}
	}
}

//load index from index file
void pyfastx_load_index(pyfastx_Index *self){
	int ret;
	sqlite3_stmt *stmt;

	PYFASTX_SQLITE_CALL(ret = sqlite3_open(self->index_file, &self->index_db));

	if (ret != SQLITE_OK) {
		PyErr_Format(PyExc_ConnectionError, "could not load index from file %s", self->index_file);
		return;
	}

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, "SELECT * FROM seq LIMIT 1;", -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	if (ret != SQLITE_ROW) {
		PyErr_Format(PyExc_RuntimeError, "the index file %s was damaged", self->index_file);
		return;
	}

	if (self->gzip_format) {
		pyfastx_load_gzip_index(self->gzip_index, self->index_db);
	}
}

void pyfastx_build_index(pyfastx_Index *self){
	PyObject *index_obj;
	index_obj = PyUnicode_FromString(self->index_file);

	if (file_exists(index_obj)) {
		pyfastx_load_index(self);
	} else {
		pyfastx_create_index(self);
	}

	Py_DECREF(index_obj);
}

void pyfastx_index_free(pyfastx_Index *self){
	if (self->gzip_format && self->gzip_index) {
		zran_free(self->gzip_index);
	}

	if (self->index_file) {
		free(self->index_file);
	}

	if (self->iter_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->iter_stmt));
	}

	if (self->uid_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->uid_stmt));
	}

	if (self->seq_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->seq_stmt));
	}

	if (self->index_db) {
		PYFASTX_SQLITE_CALL(sqlite3_close(self->index_db));
		self->index_db = NULL;
	}

	if (self->cache_seq.m) {
		free(self->cache_seq.s);
	}

	if (self->cache_name.m) {
		free(self->cache_name.s);
	}

	self->fasta = NULL;

	kseq_destroy(self->kseqs);
	fclose(self->fd);
	gzclose(self->gzfd);
}

pyfastx_Sequence* pyfastx_index_new_seq(pyfastx_Index *self) {
	pyfastx_Sequence *seq = (pyfastx_Sequence *)PyObject_CallObject((PyObject *)&pyfastx_SequenceType, NULL);
	//is full sequence
	seq->complete = 1;

	//position
	seq->start = 1;
	seq->end = seq->seq_len;

	//index
	seq->index = self;
	Py_INCREF(self->fasta);

	//buff
	seq->line_cache = NULL;
	seq->cache_pos = NULL;

	//line string init
	kstring_init(seq->line);

	seq->desc = NULL;
	seq->raw = NULL;

	return seq;
}

PyObject *pyfastx_index_make_seq(pyfastx_Index *self, sqlite3_stmt *stmt){
	//int32_t a, c, g, t, n;
	int nbytes;

	pyfastx_Sequence *seq = pyfastx_index_new_seq(self);

	PYFASTX_SQLITE_CALL(
		seq->id = sqlite3_column_int64(stmt, 0);
		nbytes = sqlite3_column_bytes(stmt, 1);
		seq->name = (char *)malloc(nbytes + 1);
		memcpy(seq->name, sqlite3_column_text(stmt, 1), nbytes);
		seq->name[nbytes] = '\0';
		seq->offset = sqlite3_column_int64(stmt, 2);
		seq->byte_len = sqlite3_column_int64(stmt, 3);
		seq->seq_len = sqlite3_column_int64(stmt, 4);
		seq->line_len = sqlite3_column_int64(stmt, 5);
		seq->end_len = sqlite3_column_int(stmt, 6);
		seq->normal = sqlite3_column_int(stmt, 7);
		seq->desc_len = sqlite3_column_int(stmt, 8);
	);

	return (PyObject *)seq;
}

PyObject *pyfastx_index_get_seq_by_name(pyfastx_Index *self, PyObject *sname){
	// sqlite3 prepare object
	int ret;
	char *name;

	Py_ssize_t nbytes;
	pyfastx_Sequence *obj;

	name = (char *)PyUnicode_AsUTF8AndSize(sname, &nbytes);

	PYFASTX_SQLITE_CALL(
		sqlite3_bind_text(self->seq_stmt, 1, name, -1, NULL);
		ret = sqlite3_step(self->seq_stmt);
	);

	if (ret == SQLITE_ROW) {
		obj = pyfastx_index_new_seq(self);
		obj->name = (char *)malloc(nbytes + 1);
		memcpy(obj->name, name, nbytes);
		obj->name[nbytes] = '\0';

		PYFASTX_SQLITE_CALL(
			obj->id = sqlite3_column_int64(self->seq_stmt, 0);
			obj->offset = sqlite3_column_int64(self->seq_stmt, 2);
			obj->byte_len = sqlite3_column_int64(self->seq_stmt, 3);
			obj->seq_len = sqlite3_column_int64(self->seq_stmt, 4);
			obj->line_len = sqlite3_column_int64(self->seq_stmt, 5);
			obj->end_len = sqlite3_column_int(self->seq_stmt, 6);
			obj->normal = sqlite3_column_int(self->seq_stmt, 7);
			obj->desc_len = sqlite3_column_int(self->seq_stmt, 8);
			sqlite3_reset(self->seq_stmt);
		);

		return (PyObject *)obj;
	} else {
		PYFASTX_SQLITE_CALL(sqlite3_reset(self->seq_stmt));
		PyErr_Format(PyExc_KeyError, "%s does not exist in fasta file", name);
		return NULL;
	}
}

PyObject *pyfastx_index_get_seq_by_id(pyfastx_Index *self, Py_ssize_t chrom){
	int ret;
	Py_ssize_t nbytes;
	pyfastx_Sequence *obj;

	PYFASTX_SQLITE_CALL(
		sqlite3_bind_int64(self->uid_stmt, 1, chrom);
		ret = sqlite3_step(self->uid_stmt);
	);

	if (ret == SQLITE_ROW){
		obj = pyfastx_index_new_seq(self);
		obj->id = chrom;
		PYFASTX_SQLITE_CALL(
			nbytes = sqlite3_column_bytes(self->uid_stmt, 1);
			obj->name = (char *)malloc(nbytes + 1);
			memcpy(obj->name, sqlite3_column_text(self->uid_stmt, 1), nbytes);
			obj->name[nbytes] = '\0';
			obj->offset = sqlite3_column_int64(self->uid_stmt, 2);
			obj->byte_len = sqlite3_column_int64(self->uid_stmt, 3);
			obj->seq_len = sqlite3_column_int64(self->uid_stmt, 4);
			obj->line_len = sqlite3_column_int64(self->uid_stmt, 5);
			obj->end_len = sqlite3_column_int(self->uid_stmt, 6);
			obj->normal = sqlite3_column_int(self->uid_stmt, 7);
			obj->desc_len = sqlite3_column_int(self->uid_stmt, 8);
			sqlite3_reset(self->uid_stmt);
		);

		return (PyObject *)obj;
	} else {
		PYFASTX_SQLITE_CALL(sqlite3_reset(self->uid_stmt));
		PyErr_SetString(PyExc_IndexError, "Index Error");
		return NULL;
	}
}

PyObject *pyfastx_index_next_null(pyfastx_Index *self) {
	PyErr_SetString(PyExc_TypeError, "'Fasta' object is not iterable");
	return NULL;
}

PyObject *pyfastx_index_next_seq(pyfastx_Index *self) {
	if (kseq_read(self->kseqs) >= 0) {
		return Py_BuildValue("(ss)", self->kseqs->name.s, self->kseqs->seq.s);
	}

	return NULL;
}

PyObject *pyfastx_index_next_upper_seq(pyfastx_Index *self) {
	if (kseq_read(self->kseqs) >= 0) {
		upper_string(self->kseqs->seq.s, self->kseqs->seq.l);
		return Py_BuildValue("(ss)", self->kseqs->name.s, self->kseqs->seq.s);
	}

	return NULL;
}

PyObject *pyfastx_index_next_full_name_seq(pyfastx_Index *self) {
	PyObject *fname;
	PyObject *ret;

	if (kseq_read(self->kseqs) >= 0) {
		if (self->kseqs->comment.l) {
			fname = PyUnicode_FromFormat("%s %s", self->kseqs->name.s, self->kseqs->comment.s);
			ret = Py_BuildValue("(Os)", fname, self->kseqs->seq.s);
			Py_DECREF(fname);
		} else {
			ret = Py_BuildValue("(ss)", self->kseqs->name.s, self->kseqs->seq.s);
		}

		return ret;
	}

	return NULL;
}

PyObject *pyfastx_index_next_full_name_upper_seq(pyfastx_Index *self) {
	PyObject *fname;
	PyObject *ret;

	if (kseq_read(self->kseqs) >= 0) {
		upper_string(self->kseqs->seq.s, self->kseqs->seq.l);

		if (self->kseqs->comment.l) {
			fname = PyUnicode_FromFormat("%s %s", self->kseqs->name.s, self->kseqs->comment.s);
			ret = Py_BuildValue("(Os)", fname, self->kseqs->seq.s);
			Py_DECREF(fname);
		} else {
			ret = Py_BuildValue("(ss)", self->kseqs->name.s, self->kseqs->seq.s);
		}

		return ret;
	}

	return NULL;
}

PyObject *pyfastx_index_next_with_index_seq(pyfastx_Index *self) {
	int ret;

	PYFASTX_SQLITE_CALL(ret = sqlite3_step(self->iter_stmt));

	if (ret == SQLITE_ROW) {
		return pyfastx_index_make_seq(self, self->iter_stmt);
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(self->iter_stmt));

	self->iterating = 0;
	self->iter_stmt = NULL;

	return NULL;
}

void pyfastx_index_random_read(pyfastx_Index* self, char* buff, Py_ssize_t offset, Py_ssize_t bytes) {
	if (self->gzip_format) {
		zran_seek(self->gzip_index, offset, SEEK_SET, NULL);
		zran_read(self->gzip_index, buff, bytes);
	} else {
		FSEEK(self->fd, offset, SEEK_SET);
		fread(buff, bytes, 1, self->fd);
	}
	buff[bytes] = '\0';
}

void pyfastx_index_fill_cache(pyfastx_Index* self, Py_ssize_t offset, Py_ssize_t size) {
	if (size >= self->cache_seq.m ) {
		self->cache_seq.m = size + 1;
		self->cache_seq.s = (char *)realloc(self->cache_seq.s, self->cache_seq.m);
	}

	pyfastx_index_random_read(self, self->cache_seq.s, offset, size);

	if (self->uppercase) {
		self->cache_seq.l = remove_space_uppercase(self->cache_seq.s, size);
	} else {
		self->cache_seq.l = remove_space(self->cache_seq.s, size);
	}
}
