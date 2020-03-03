#include "index.h"
#include "util.h"
#include "kseq.h"
#include "sequence.h"

/*
create a index
@param file_path, fasta path and name
@param uppercase, uppercase sequence
@param uppercase
*/
pyfastx_Index* pyfastx_init_index(char* file_name, int file_len, int uppercase, int memory_index, PyObject* key_func){
	pyfastx_Index* index;

	index = (pyfastx_Index *)malloc(sizeof(pyfastx_Index));
	index->uppercase = uppercase;

	//key function
	index->key_func = key_func;

	//check input file is gzip or not
	index->gzip_format = is_gzip_format(file_name);

	//initial kseqs
	index->gzfd = gzopen(file_name, "rb");
	index->kseqs = kseq_init(index->gzfd);

	//create index file or memory index
	if (memory_index) {
		index->index_file = ":memory:";
	} else {
		index->index_file = (char *)malloc(file_len + 5);
		strcpy(index->index_file, file_name);
		strcat(index->index_file, ".fxi");
	}

	index->fd = fopen(file_name, "rb");

	index->index_db = 0;

	if(index->gzip_format){
		index->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
		//initial zran index
		zran_init(index->gzip_index, index->fd, 4194304, 32768, 1048576, ZRAN_AUTO_BUILD);
	}

	//cache name
	index->cache_chrom = 0;

	//cache start and end position
	index->cache_start = 0;
	index->cache_end = 0;

	//cache sequence
	index->cache_seq = NULL;

	return index;
}

void pyfastx_rewind_index(pyfastx_Index *self){
	kseq_rewind(self->kseqs);
	gzrewind(self->gzfd);
}

PyObject* pyfastx_get_next_seq(pyfastx_Index *self){
	if (kseq_read(self->kseqs) >= 0) {
		if (self->uppercase) {
			upper_string(self->kseqs->seq.s);
		}
		return Py_BuildValue("(ss)", self->kseqs->name.s, self->kseqs->seq.s);
	}
	return NULL;
}

void pyfastx_create_index(pyfastx_Index *self){
	// seqlite3 return value
	int ret;
	
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	
	// 1: normal fasta sequence with the same length in line
	// 0: not normal fasta sequence with different length in line
	uint16_t seq_normal = 1;

	//1: \n, 2: \r\n
	uint16_t line_end = 1;
	
	//length of previous line
	uint32_t line_len = 0;

	//length of current line
	uint32_t temp_len = 0;

	//number of lines that line_len not equal to temp_len
	uint32_t bad_line = 0;

	//current read position
	uint64_t position = 0;
	
	// start position
	uint64_t start = 0;

	//sequence length
	uint32_t seq_len = 0;

	//real line len
	uint32_t real_len;

	//total sequence count
	uint32_t total_seq = 0;

	//total sequence length
	uint64_t total_len = 0;

	//space position
	char *space_pos;

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
		PyErr_Format(PyExc_ConnectionError, "Can not open index file %s", self->index_file);
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
			llen INTEGER, --line lenght\n \
			elen INTEGER, --end length\n \
			norm INTEGER, --line with the same length or not\n \
			dlen INTEGER --description header line length\n \
		); \
		CREATE TABLE stat ( \
			seqnum INTEGER, --total seq counts \n \
			seqlen INTEGER --total seq length \n \
		); \
		CREATE TABLE comp ( \
			ID INTEGER PRIMARY KEY, \
			a INTEGER, \
			b INTEGER, \
			c INTEGER, \
			d INTEGER, \
			e INTEGER, \
			f INTEGER, \
			g INTEGER, \
			h INTEGER, \
			i INTEGER, \
			j INTEGER, \
			k INTEGER, \
			l INTEGER, \
			m INTEGER, \
			n INTEGER, \
			o INTEGER, \
			p INTEGER, \
			q INTEGER, \
			r INTEGER, \
			s INTEGER, \
			t INTEGER, \
			u INTEGER, \
			v INTEGER, \
			w INTEGER, \
			x INTEGER, \
			y INTEGER, \
			z INTEGER \
		); \
		CREATE TABLE gzindex ( \
			ID INTEGER PRIMARY KEY, \
			content BLOB \
		);";

	PYFASTX_SQLITE_CALL(ret = sqlite3_exec(self->index_db, sql, NULL, NULL, NULL));

	if (ret != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, "Can not create index tables");
		return;
	}
	
	sql = "PRAGMA synchronous=OFF;BEGIN TRANSACTION;";
	PYFASTX_SQLITE_CALL(ret = sqlite3_exec(self->index_db, sql, NULL, NULL, NULL));
	if (ret != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, "Can not begin transaction");
		return;
	}

	sql = "INSERT INTO seq VALUES (?,?,?,?,?,?,?,?,?);";
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));
	
	gzrewind(self->gzfd);
	ks = ks_init(self->gzfd);

	Py_BEGIN_ALLOW_THREADS
	while (ks_getuntil(ks, '\n', &line, 0) >= 0) {
		position += line.l + 1;

		//first char is >
		if (line.s[0] == 62) {
			if (start > 0) {
				//end of sequence and check whether normal fasta
				seq_normal = (bad_line > 1) ? 0 : 1;
				
				sqlite3_bind_null(stmt, 1);
				sqlite3_bind_text(stmt, 2, chrom.s, chrom.l, SQLITE_STATIC);
				sqlite3_bind_int64(stmt, 3, start);
				sqlite3_bind_int(stmt, 4, position-start-line.l-1);
				sqlite3_bind_int(stmt, 5, seq_len);
				sqlite3_bind_int(stmt, 6, line_len);
				sqlite3_bind_int(stmt, 7, line_end);
				sqlite3_bind_int(stmt, 8, seq_normal);
				sqlite3_bind_int(stmt, 9, desc_len);
				sqlite3_step(stmt);
				sqlite3_reset(stmt);

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

			if (self->key_func == Py_None) {
				space_pos = strchr(header_pos, ' ');
				if (space_pos == NULL) {
					chrom.l = desc_len;
				} else {
					chrom.l = space_pos - header_pos;
				}
				memcpy(chrom.s, header_pos, chrom.l);
				chrom.s[chrom.l] = '\0';
			} else {
				PyGILState_STATE state = PyGILState_Ensure();
				PyObject *result = PyObject_CallFunction(self->key_func, "s", header_pos);
				PyGILState_Release(state);
				temp_chrom = PyUnicode_AsUTF8AndSize(result, &chrom.l);
				memcpy(chrom.s, temp_chrom, chrom.l);
				chrom.s[chrom.l] = '\0';
				free(temp_chrom);
				Py_DECREF(result);
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
	
	sqlite3_bind_null(stmt, 1);
	sqlite3_bind_text(stmt, 2, chrom.s, chrom.l, SQLITE_STATIC);
	sqlite3_bind_int64(stmt, 3, start);
	sqlite3_bind_int(stmt, 4, position-start);
	sqlite3_bind_int(stmt, 5, seq_len);
	sqlite3_bind_int(stmt, 6, line_len);
	sqlite3_bind_int(stmt, 7, line_end);
	sqlite3_bind_int(stmt, 8, seq_normal);
	sqlite3_bind_int(stmt, 9, desc_len);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	stmt = NULL;

	++total_seq;
	total_len += seq_len;
	
	sqlite3_exec(self->index_db, "COMMIT;", NULL, NULL, NULL);
	sqlite3_exec(self->index_db, "CREATE INDEX chromidx ON seq (chrom);", NULL, NULL, NULL);

	sqlite3_prepare_v2(self->index_db, "INSERT INTO stat VALUES (?,?);", -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, total_seq);
	sqlite3_bind_int64(stmt, 2, total_len);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	ks_destroy(ks);
	free(line.s);
	free(chrom.s);

	//create gzip random access index
	if (self->gzip_format) {
		if (strcmp(self->index_file, ":memory:") == 0) {
			zran_build_index(self->gzip_index, 0, 0);
		} else {
			pyfastx_build_gzip_index(self->gzip_index, self->index_db, self->index_file);
		}
	}

	Py_END_ALLOW_THREADS
}

//load index from index file
void pyfastx_load_index(pyfastx_Index *self){
	int ret;
	sqlite3_stmt *stmt;

	PYFASTX_SQLITE_CALL(ret = sqlite3_open(self->index_file, &self->index_db));

	if (ret != SQLITE_OK) {
		PyErr_Format(PyExc_ConnectionError, "can not load index from file %s", self->index_file);
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
		pyfastx_load_gzip_index(self->gzip_index, self->index_db, self->index_file);
	}
}

void pyfastx_build_index(pyfastx_Index *self){
	if (file_exists(self->index_file)) {
		pyfastx_load_index(self);
	} else {
		pyfastx_create_index(self);
	}
}

void pyfastx_index_free(pyfastx_Index *self){
	if (self->gzip_format && self->gzip_index != NULL) {
		zran_free(self->gzip_index);
	}

	if (self->index_db) {
		PYFASTX_SQLITE_CALL(sqlite3_close(self->index_db));
	}

	if (self->cache_seq) {
		free(self->cache_seq);
	}

	kseq_destroy(self->kseqs);
	fclose(self->fd);
	gzclose(self->gzfd);
}

PyObject *pyfastx_index_make_seq(pyfastx_Index *self, sqlite3_stmt *stmt){
	//int32_t a, c, g, t, n;
	int nbytes;

	//pyfastx_Sequence *seq = PyObject_New(pyfastx_Sequence, &pyfastx_SequenceType);
	pyfastx_Sequence *seq = (pyfastx_Sequence *)PyObject_CallObject((PyObject *)&pyfastx_SequenceType, NULL);

	if(!seq){
		return NULL;
	}

	PYFASTX_SQLITE_CALL(
		seq->id = sqlite3_column_int(stmt, 0);
		nbytes = sqlite3_column_bytes(stmt, 1);
		seq->name = (char *)malloc(nbytes + 1);
		memcpy(seq->name, (char *)sqlite3_column_text(stmt, 1), nbytes);
		seq->name[nbytes] = '\0';
		seq->offset = (int64_t)sqlite3_column_int64(stmt, 2);
		seq->byte_len = sqlite3_column_int(stmt, 3);
		seq->seq_len = sqlite3_column_int(stmt, 4);
		seq->line_len = sqlite3_column_int(stmt, 5);
		seq->end_len = sqlite3_column_int(stmt, 6);
		seq->normal = sqlite3_column_int(stmt, 7);
		//sqlite3_finalize(stmt);
	);
	
	seq->parent_len = seq->seq_len;
	seq->ks = NULL;

	//position
	seq->start = 1;
	seq->end = seq->seq_len;

	//index
	seq->index = self;

	//Py_INCREF(seq);
	return (PyObject *)seq;
}

PyObject *pyfastx_index_get_seq_by_name(pyfastx_Index *self, char *name){
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	PyObject *obj;
	int ret;
	
	//select sql statement, chrom indicates seq name or chromomsome
	const char* sql = "SELECT * FROM seq WHERE chrom=? LIMIT 1;";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
	);
	
	if(ret != SQLITE_ROW){
		PyErr_Format(PyExc_KeyError, "%s does not exist in fasta file", name);
		return NULL;
	}

	obj = pyfastx_index_make_seq(self, stmt);
	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return obj;
}


PyObject *pyfastx_index_get_seq_by_id(pyfastx_Index *self, uint32_t chrom){
	sqlite3_stmt *stmt;
	PyObject *obj;
	int ret;

	const char* sql = "SELECT * FROM seq WHERE ID=? LIMIT 1;";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int(stmt, 1, chrom);
		ret = sqlite3_step(stmt);
	);
	
	if (ret != SQLITE_ROW){
		PyErr_SetString(PyExc_IndexError, "Index Error");
		return NULL;
	}

	obj = pyfastx_index_make_seq(self, stmt);
	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return obj;
}

char *pyfastx_index_get_full_seq(pyfastx_Index *self, uint32_t chrom){
	sqlite3_stmt *stmt;
	uint32_t seq_len;
	int64_t offset;
	uint32_t bytes;
	int ret;
	//char *buff;

	//select sql statement, chrom indicates seq or chromomsome id
	const char* sql = "SELECT boff,blen,slen FROM seq WHERE ID=? LIMIT 1;";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int(stmt, 1, chrom);
		ret = sqlite3_step(stmt);
	);

	if(ret != SQLITE_ROW){
		PyErr_SetString(PyExc_KeyError, "Can not found sequence");
		return NULL;
	}

	PYFASTX_SQLITE_CALL(
		offset = (int64_t)sqlite3_column_int64(stmt, 0);
		bytes = sqlite3_column_int(stmt, 1);
		seq_len = sqlite3_column_int(stmt, 2);
		sqlite3_finalize(stmt);
	);

	if((chrom == self->cache_chrom) && (1==self->cache_start) && (seq_len==self->cache_end)){
		return self->cache_seq;
	}

	//Py_BEGIN_ALLOW_THREADS

	self->cache_seq = (char *)malloc(bytes + 1);
	
	if (self->gzip_format) {
		zran_seek(self->gzip_index, offset, SEEK_SET, NULL);
		zran_read(self->gzip_index, self->cache_seq, bytes);
	} else {
		fseek(self->fd, offset, SEEK_SET);
		if (fread(self->cache_seq, bytes, 1, self->fd) != 1) {
			if (!feof(self->fd)) {
				PyErr_SetString(PyExc_RuntimeError, "reading sequence error");
				return NULL;
			}
		}
	}

	self->cache_seq[bytes] = '\0';

	if (self->uppercase) {
		remove_space_uppercase(self->cache_seq);
	} else {
		remove_space(self->cache_seq);
	}

	//Py_END_ALLOW_THREADS

	self->cache_chrom = chrom;
	self->cache_start = 1;
	self->cache_end = seq_len;

	return self->cache_seq;
}

//clear index cached sequence
void pyfastx_index_cache_clear(pyfastx_Index *self) {
	self->cache_chrom = 0;
	self->cache_start = 0;
	self->cache_end = 0;
	self->cache_seq = NULL;
}
