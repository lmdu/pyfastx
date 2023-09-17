#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "fastq.h"
#include "read.h"
#include "fqkeys.h"
#include "structmember.h"

void pyfastx_fastq_create_index(pyfastx_Fastq *self) {
	int ret;
	int l;
	int j;
	int dlen = 0;

	char* space;
	const char *sql;

	sqlite3_stmt *stmt;

	Py_ssize_t rlen = 0;
	Py_ssize_t soff = 0;
	Py_ssize_t qoff = 0;
	Py_ssize_t pos = 0;
	Py_ssize_t size = 0;
	Py_ssize_t line_num = 0;
	
	kstring_t name = {0, 0, 0};
	kstring_t line = {0, 0, 0};

	sql = " \
		CREATE TABLE read ( \
			ID INTEGER PRIMARY KEY, --read id \n \
			name TEXT, --read name \n \
			dlen INTEGER, --description length \n \
			rlen INTEGER, --read length \n \
			soff INTEGER, --read seq offset \n \
			qoff INTEGER --read qual offset \n \
		); \
		CREATE TABLE gzindex (  \
			ID INTEGER PRIMARY KEY,  \
			content BLOB  \
		); \
		CREATE TABLE stat ( \
			counts INTEGER, --read counts \n \
			size INTEGER, --all read length \n \
			avglen REAL --average read length \n \
		); \
		CREATE TABLE base ( \
			a INTEGER,  \
			c INTEGER,  \
			g INTEGER,  \
			t INTEGER,  \
			n INTEGER  \
		); \
		CREATE TABLE meta ( \
			maxlen INTEGER, --maximum read length \n \
			minlen INTEGER, --minimum read length \n \
			minqs INTEGER, --max quality score \n \
			maxqs INTEGER, --min quality score \n \
			phred INTEGER --phred value \n \
		);";

	PYFASTX_SQLITE_CALL(ret=sqlite3_open(self->index_file, &self->index_db));
	if (ret != SQLITE_OK){
		PyErr_Format(PyExc_ConnectionError, "could not open index file %s", self->index_file);
		return;
	}

	PYFASTX_SQLITE_CALL(ret=sqlite3_exec(self->index_db, sql, NULL, NULL, NULL));
	if (ret != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, "could not create index table");
		return;
	}

	sql = "PRAGMA synchronous = OFF; PRAGMA locking_mode=EXCLUSIVE; BEGIN TRANSACTION;";
	PYFASTX_SQLITE_CALL(ret=sqlite3_exec(self->index_db, sql, NULL, NULL, NULL));
	if(ret != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, "can not begin transaction");
		return;
	}
	
	sql = "INSERT INTO read VALUES (?,?,?,?,?,?);";
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));

	gzrewind(self->middle->gzfd);
	ks_rewind(self->ks);

	//Py_BEGIN_ALLOW_THREADS

	while ((l=ks_getuntil(self->ks, '\n', &line, 0)) >= 0) {
		++line_num;
		++l;

		j = line_num % 4;

		switch(j) {
			case 1:
				//name.m max length of name line in fastq file
				//name.l real length of name
				if (name.m < line.l) {
					name.s = (char *)realloc(name.s, line.l);
					name.m = line.l;
				}
				dlen = line.l;
				name.l = line.l-1;
				memcpy(name.s, line.s+1, name.l);

				if (name.s[name.l-1] == '\r') {
					--name.l;
				}

				space = strchr(name.s, ' ');

				if (space != NULL) {
					name.l = space - name.s;
				}

				break;

			case 2:
				soff = pos;
				
				if (line.s[line.l-1] == '\r') {
					rlen = line.l - 1;
				} else {
					rlen = line.l;
				}
				size += rlen;
				break;

			case 0:
				qoff = pos;

				//write to sqlite3
				PYFASTX_SQLITE_CALL(
					sqlite3_bind_null(stmt, 1);
					sqlite3_bind_text(stmt, 2, name.s, name.l, SQLITE_STATIC);
					sqlite3_bind_int(stmt, 3, dlen);
					sqlite3_bind_int64(stmt, 4, rlen);
					sqlite3_bind_int64(stmt, 5, soff);
					sqlite3_bind_int64(stmt, 6, qoff);
					sqlite3_step(stmt);
					sqlite3_reset(stmt);
				);
				break;
		}
		pos += l;
	}

	PYFASTX_SQLITE_CALL(
		sqlite3_finalize(stmt);
		sqlite3_exec(self->index_db, "PRAGMA locking_mode=NORMAL;", NULL, NULL, NULL);
		sqlite3_exec(self->index_db, "COMMIT;", NULL, NULL, NULL);
		sqlite3_exec(self->index_db, "CREATE UNIQUE INDEX readidx ON read (name);", NULL, NULL, NULL);
	);
	stmt = NULL;

	self->read_counts = line_num/4;
	self->seq_length = size;
	self->avg_length = size*1.0/self->read_counts;
	sql = "INSERT INTO stat VALUES (?,?,?);";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int64(stmt, 1, self->read_counts);
		sqlite3_bind_int64(stmt, 2, self->seq_length);
		sqlite3_bind_double(stmt, 3, self->avg_length);
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	//Py_END_ALLOW_THREADS

	free(line.s);
	free(name.s);

	//if is gzip build gzip index
	if (self->middle->gzip_format) {
		pyfastx_build_gzip_index(self->middle->gzip_index, self->index_db);
	}
}

void pyfastx_fastq_load_index(pyfastx_Fastq *self) {
	int ret;
	const char* sql;
	sqlite3_stmt* stmt;	

	PYFASTX_SQLITE_CALL(ret=sqlite3_open(self->index_file, &self->index_db));

	if (ret != SQLITE_OK){
		PyErr_Format(PyExc_ConnectionError, "could not open index file %s", self->index_file);
		return;
	}

	//calculate attributes
	sql = "SELECT * FROM stat LIMIT 1;";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			self->read_counts = sqlite3_column_int64(stmt, 0);
			self->seq_length = sqlite3_column_int64(stmt, 1);
			self->avg_length = sqlite3_column_double(stmt, 2);
			sqlite3_finalize(stmt);
		);
	} else {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_Format(PyExc_RuntimeError, "the index file %s was damaged", self->index_file);
		return;
	}
	
	stmt = NULL;

	sql = "SELECT phred FROM meta LIMIT 1;";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(self->middle->phred = sqlite3_column_int(stmt, 0));
	} else {
		self->middle->phred = 0;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	if(self->middle->gzip_format){
		pyfastx_load_gzip_index(self->middle->gzip_index, self->index_db);
	}
}

PyObject *pyfastx_fastq_build_index(pyfastx_Fastq *self) {
	PyObject *index_obj;
	index_obj = PyUnicode_FromString(self->index_file);

	if (file_exists(index_obj)) {
		pyfastx_fastq_load_index(self);
	} else {
		pyfastx_fastq_create_index(self);
	}

	Py_DECREF(index_obj);

	Py_RETURN_TRUE;
}

PyObject *pyfastx_fastq_next_null(pyfastx_FastqMiddleware *self) {
	PyErr_SetString(PyExc_TypeError, "'Fastq' object is not iterable");
	return NULL;
}

PyObject *pyfastx_fastq_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	int phred = 0;
	int build_index = 1;
	int full_index = 0;
	int full_name = 0;

	char *index_file;

	PyObject *file_obj;
	PyObject *index_obj = NULL;

	Py_ssize_t index_len;

	static char* keywords[] = {"file_name", "index_file", "phred", "build_index", "full_index", "full_name", NULL};

	pyfastx_Fastq *obj;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|Oiiii", keywords, &file_obj, &index_obj, &phred, &build_index, &full_index, &full_name)) {
		return NULL;
	}

	if (!file_exists(file_obj)) {
		PyErr_Format(PyExc_FileExistsError, "input fastq file %U does not exists", file_obj);
		return NULL;
	}

	obj = (pyfastx_Fastq *)type->tp_alloc(type, 0);
	if (!obj) {
		return NULL;
	}

	obj->middle = (pyfastx_FastqMiddleware *)malloc(sizeof(pyfastx_FastqMiddleware));

	obj->file_obj = Py_NewRef(file_obj);

	//check input file is gzip or not
	obj->middle->gzip_format = is_gzip_format(file_obj);

	//initial kstream and kseq
	obj->middle->gzfd = pyfastx_gzip_open(file_obj, "rb");
	obj->ks = ks_init(obj->middle->gzfd);
	obj->middle->kseq = kseq_init(obj->middle->gzfd);

	//check is correct fastq format
	if (!fastq_validator(obj->middle->gzfd)) {
		PyErr_Format(PyExc_RuntimeError, "%U is not plain or gzip compressed fastq formatted file", file_obj);
		return NULL;
	}

	//create index file
	if (index_obj) {
		index_file = (char *)PyUnicode_AsUTF8AndSize(index_obj, &index_len);
		obj->index_file = (char *)malloc(index_len);
		memcpy(obj->index_file, index_file, index_len);
		obj->index_file[index_len] = '\0';
	} else {
		index_file = (char *)PyUnicode_AsUTF8AndSize(file_obj, &index_len);
		index_len += 5;
		obj->index_file = (char *)malloc(index_len);
		strcpy(obj->index_file, index_file);
		strcat(obj->index_file, ".fxi");
	}

	obj->middle->fd = _Py_fopen_obj(obj->file_obj, "rb");

	//initail index connection
	obj->index_db = 0;
	obj->middle->iter_stmt = NULL;
	obj->id_stmt = NULL;
	obj->name_stmt = NULL;

	obj->has_index = build_index;
	obj->full_name = full_name;

	//initialize attribute
	obj->middle->phred = phred;
	obj->gc_content = 0;
	obj->minlen = 0;
	obj->maxlen = 0;
	obj->minqual = 0;
	obj->maxqual = 0;

	//is gzip file
	if (obj->middle->gzip_format) {
		obj->middle->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
		zran_init(obj->middle->gzip_index, obj->middle->fd, NULL, 1048576, 32768, 16384, ZRAN_AUTO_BUILD);
	}

	index_obj = PyUnicode_FromString(obj->index_file);

	if (file_exists(index_obj)) {
		pyfastx_fastq_load_index(obj);
	} else if (build_index) {
		pyfastx_fastq_create_index(obj);
	}

	Py_DECREF(index_obj);

	//prepare sql
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(obj->index_db, "SELECT * FROM read WHERE ID=? LIMIT 1", -1, &obj->id_stmt, NULL);
		sqlite3_prepare_v2(obj->index_db, "SELECT * FROM read WHERE name=? LIMIT 1", -1, &obj->name_stmt, NULL);
	);

	if (build_index && full_index) {
		pyfastx_fastq_calc_composition(obj);
	}

	//iter function
	obj->func = pyfastx_fastq_next_null;

	//initialize cache buffer
	obj->middle->cache_buff = NULL;
	obj->middle->cache_soff = 0;
	obj->middle->cache_eoff = 0;

	obj->middle->fastq = (PyObject *)obj;

	return (PyObject *)obj;
}

void pyfastx_fastq_dealloc(pyfastx_Fastq *self) {
	if (self->middle->iter_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->middle->iter_stmt));
	}

	if (self->id_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->id_stmt));
	}

	if (self->name_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->name_stmt));
	}

	if (self->index_db) {
		PYFASTX_SQLITE_CALL(sqlite3_close(self->index_db));
	}

	if (self->middle->gzip_format) {
		zran_free(self->middle->gzip_index);
	}

	if (self->middle->cache_buff) {
		free(self->middle->cache_buff);
	}

	self->middle->fastq = NULL;

	ks_destroy(self->ks);
	kseq_destroy(self->middle->kseq);
	fclose(self->middle->fd);
	gzclose(self->middle->gzfd);

	Py_DECREF(self->file_obj);

	Py_TYPE(self)->tp_free((PyObject *)self);
}

Py_ssize_t pyfastx_fastq_length(pyfastx_Fastq *self) {
	return self->read_counts;
}

pyfastx_Read* pyfastx_fastq_new_read(pyfastx_FastqMiddleware *middle) {
	pyfastx_Read* read = (pyfastx_Read *)PyObject_CallObject((PyObject *)&pyfastx_ReadType, NULL);
	read->middle = middle;
	read->seq = NULL;
	read->qual = NULL;
	read->raw = NULL;
	read->desc = NULL;

	Py_INCREF(middle->fastq);

	return read;
}

PyObject* pyfastx_fastq_make_read(pyfastx_FastqMiddleware *middle) {
	int nbytes;

	//pyfastx_Read *read = PyObject_New(pyfastx_Read, &pyfastx_ReadType);
	pyfastx_Read *read = pyfastx_fastq_new_read(middle);

	PYFASTX_SQLITE_CALL(
		read->id = sqlite3_column_int64(middle->iter_stmt, 0);
		nbytes = sqlite3_column_bytes(middle->iter_stmt, 1);
		read->name = (char *)malloc(nbytes + 1);
		memcpy(read->name, sqlite3_column_text(middle->iter_stmt, 1), nbytes);
		read->name[nbytes] = '\0';
		read->desc_len = sqlite3_column_int(middle->iter_stmt, 2);
		read->read_len = sqlite3_column_int64(middle->iter_stmt, 3);
		read->seq_offset = sqlite3_column_int64(middle->iter_stmt, 4);
		read->qual_offset = sqlite3_column_int64(middle->iter_stmt, 5);
		//sqlite3_finalize(stmt);
	);

	return (PyObject *)read;
}

PyObject* pyfastx_fastq_get_read_by_id(pyfastx_Fastq *self, Py_ssize_t read_id) {
	int ret;
	int nbytes;
	pyfastx_Read *obj;

	PYFASTX_SQLITE_CALL(
		sqlite3_bind_int(self->id_stmt, 1, read_id);
		ret = sqlite3_step(self->id_stmt);
	);

	if (ret == SQLITE_ROW) {
		obj = pyfastx_fastq_new_read(self->middle);
		obj->id = read_id;
		PYFASTX_SQLITE_CALL(
			nbytes = sqlite3_column_bytes(self->id_stmt, 1);
			obj->name = (char *)malloc(nbytes + 1);
			memcpy(obj->name, sqlite3_column_text(self->id_stmt, 1), nbytes);
			obj->name[nbytes] = '\0';
			obj->desc_len = sqlite3_column_int(self->id_stmt, 2);
			obj->read_len = sqlite3_column_int64(self->id_stmt, 3);
			obj->seq_offset = sqlite3_column_int64(self->id_stmt, 4);
			obj->qual_offset = sqlite3_column_int64(self->id_stmt, 5);
			sqlite3_reset(self->id_stmt);
		);

		return (PyObject *)obj;
	} else {
		PyErr_SetString(PyExc_IndexError, "Index Error");
		return NULL;
	}
}

PyObject* pyfastx_fastq_get_read_by_name(pyfastx_Fastq *self, PyObject* rname) {
	int ret;
	char *name;
	Py_ssize_t nbytes;
	pyfastx_Read *obj;

	name = (char *)PyUnicode_AsUTF8AndSize(rname, &nbytes);

	PYFASTX_SQLITE_CALL(
		sqlite3_bind_text(self->name_stmt, 1, name, -1, NULL);
		ret = sqlite3_step(self->name_stmt);
	);

	if (ret == SQLITE_ROW) {
		obj = pyfastx_fastq_new_read(self->middle);
		obj->name = (char *)malloc(nbytes + 1);
		memcpy(obj->name, name, nbytes);
		obj->name[nbytes] = '\0';

		PYFASTX_SQLITE_CALL(
			obj->id = sqlite3_column_int64(self->name_stmt, 0);
			obj->desc_len = sqlite3_column_int(self->name_stmt, 2);
			obj->read_len = sqlite3_column_int64(self->name_stmt, 3);
			obj->seq_offset = sqlite3_column_int64(self->name_stmt, 4);
			obj->qual_offset = sqlite3_column_int64(self->name_stmt, 5);
			sqlite3_reset(self->name_stmt);
		);

		return (PyObject *)obj;
	} else {
		PyErr_Format(PyExc_KeyError, "%s does not exist in fastq file", name);
		return NULL;
	}
}

PyObject* pyfastx_fastq_subscript(pyfastx_Fastq *self, PyObject *item) {
	Py_ssize_t i;

	self->middle->iterating = 0;

	if (PyUnicode_Check(item)) {
		return pyfastx_fastq_get_read_by_name(self, item);
	} else if (PyIndex_Check(item)) {
		i = PyNumber_AsSsize_t(item, PyExc_IndexError);

		if (i < 0) {
			i += self->read_counts;
		}

		if (i >= self->read_counts) {
			PyErr_SetString(PyExc_IndexError, "index out of range");
			return NULL;
		}

		return pyfastx_fastq_get_read_by_id(self, i+1);
	} else {
		PyErr_SetString(PyExc_KeyError, "the key must be index number or read name");
		return NULL;
	}
}

PyObject* pyfastx_fastq_repr(pyfastx_Fastq *self) {
	if (self->has_index) {
		return PyUnicode_FromFormat("<Fastq> %U contains %zd reads", self->file_obj, self->read_counts);
	} else {
		return PyUnicode_FromFormat("<Fastq> %U", self->file_obj);
	}
}

int pyfastx_fastq_contains(pyfastx_Fastq *self, PyObject *key) {
	int ret;
	char *name;
	const char *sql;
	sqlite3_stmt *stmt;

	if (!PyUnicode_Check(key)) {
		return 0;
	}

	name = (char *)PyUnicode_AsUTF8(key);
	sql = "SELECT 1 FROM read WHERE name=? LIMIT 1;";

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	return ret==SQLITE_ROW ? 1 : 0;
}

PyObject *pyfastx_fastq_next_with_index_read(pyfastx_FastqMiddleware *middle) {
	int ret;

	PYFASTX_SQLITE_CALL(ret = sqlite3_step(middle->iter_stmt));

	if (ret == SQLITE_ROW) {
		return pyfastx_fastq_make_read(middle);
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(middle->iter_stmt));
	middle->iter_stmt = NULL;
	middle->iterating = 0;

	if (middle->cache_buff) {
		free(middle->cache_buff);
		middle->cache_buff = NULL;
	}

	return NULL;
}

PyObject *pyfastx_fastq_next_read(pyfastx_FastqMiddleware *middle) {
	if (kseq_read(middle->kseq) >= 0) {
		return Py_BuildValue("sss", middle->kseq->name.s, middle->kseq->seq.s, middle->kseq->qual.s);
	}

	return NULL;
}

PyObject *pyfastx_fastq_next_full_name_read(pyfastx_FastqMiddleware *middle) {
	PyObject *rname;
	PyObject *ret;

	if (kseq_read(middle->kseq) >= 0) {
		if (middle->kseq->comment.l) {
			rname = PyUnicode_FromFormat("%s %s", middle->kseq->name.s, middle->kseq->comment.s);
			ret = Py_BuildValue("Oss", rname, middle->kseq->seq.s, middle->kseq->qual.s);
			Py_DECREF(rname);
		} else {
			ret = Py_BuildValue("sss", middle->kseq->name.s, middle->kseq->seq.s, middle->kseq->qual.s);
		}
		return ret;
	}

	return NULL;
}

PyObject *pyfastx_fastq_iter(pyfastx_Fastq *self) {
	gzrewind(self->middle->gzfd);
	rewind(self->middle->fd);
	
	if (self->has_index) {
		self->middle->iterating = 1;

		if (!self->middle->cache_buff) {
			self->middle->cache_buff = (char *)malloc(1048576);
		}
		self->middle->cache_soff = 0;
		self->middle->cache_eoff = 0;

		PYFASTX_SQLITE_CALL(
			sqlite3_finalize(self->middle->iter_stmt);
			self->middle->iter_stmt = NULL;
			sqlite3_prepare_v2(self->index_db, "SELECT * FROM read", -1, &self->middle->iter_stmt, NULL);
		);

		self->func = pyfastx_fastq_next_with_index_read;
	} else {
		kseq_rewind(self->middle->kseq);

		if (self->full_name) {
			self->func = pyfastx_fastq_next_full_name_read;
		} else {
			self->func = pyfastx_fastq_next_read;
		}
	}
	
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fastq_next(pyfastx_Fastq *self) {
	return self->func(self->middle);
}

void pyfastx_fastq_calc_composition(pyfastx_Fastq *self) {
	int i, j, ret;

	//min and max quality score
	int minqs = 104;
	int maxqs = 33;
	int phred = 0;

	const char *sql;

	sqlite3_stmt *stmt;
	kstring_t line = {0, 0, 0};

	//min and max read length
	Py_ssize_t maxlen = 0;
	Py_ssize_t minlen = 10000000000;

	//base number
	Py_ssize_t a = 0, c = 0, g = 0, t = 0, n = 0;
	Py_ssize_t line_num = 0;
	

	sql = "SELECT * FROM meta LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			if (!self->maxlen)
				self->maxlen = sqlite3_column_int64(stmt, 0);
			if (!self->minlen)
				self->minlen = sqlite3_column_int64(stmt, 1);
			if (!self->minqual)
				self->minqual = sqlite3_column_int(stmt, 2);
			if (!self->maxqual)
				self->maxqual = sqlite3_column_int(stmt, 3);
			if (!self->middle->phred)
				self->middle->phred = sqlite3_column_int(stmt, 4);
		);

		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		return;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
	stmt = NULL;

	gzrewind(self->middle->gzfd);
	ks_rewind(self->ks);

	while (ks_getuntil(self->ks, '\n', &line, 0) >= 0) {
		++line_num;

		j = line_num % 4;

		if (j == 2) {
			for (i = 0; i < line.l; i++) {
				switch (line.s[i]) {
					case 65: ++a; break;
					case 67: ++c; break;
					case 71: ++g; break;
					case 84: ++t; break;
					case 13: break;
					default: ++n;
				}
			}
		} else if (j == 0) {
			for (i = 0; i < line.l; i++) {
				if (line.s[i] == 13) {
					--line.l;
					continue;
				}

				if (line.s[i] < minqs) {
					minqs = line.s[i];
				}

				if (line.s[i] > maxqs) {
					maxqs = line.s[i];
				}
			}

			if (line.l > maxlen)
				maxlen = line.l;

			if (line.l < minlen)
				minlen = line.l;
		}
	}

	sql = "INSERT INTO base VALUES (?,?,?,?,?);";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int64(stmt, 1, a);
		sqlite3_bind_int64(stmt, 2, c);
		sqlite3_bind_int64(stmt, 3, g);
		sqlite3_bind_int64(stmt, 4, t);
		sqlite3_bind_int64(stmt, 5, n);
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);
	stmt = NULL;

	if (maxqs > 74) {
		phred = 64;
	}

	if (minqs < 59) {
		phred = 33;
	}

	//insert platform into index file
	sql = "INSERT INTO meta VALUES (?,?,?,?,?);";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int64(stmt, 1, maxlen);
		sqlite3_bind_int64(stmt, 2, minlen);
		sqlite3_bind_int(stmt, 3, minqs);
		sqlite3_bind_int(stmt, 4, maxqs);
		sqlite3_bind_int(stmt, 5, phred);
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);
	stmt = NULL;

	self->minlen = minlen;
	self->maxlen = maxlen;
	self->minqual = minqs;
	self->maxqual = maxqs;
	self->middle->phred = phred;
}

PyObject* pyfastx_fastq_guess_encoding_type(pyfastx_Fastq* self, void* closure) {
	int ret;
	int maxqs;
	int minqs;

	const char *sql;

	sqlite3_stmt *stmt;

	PyObject* platforms;
	PyObject* platform;

	pyfastx_fastq_calc_composition(self);

	sql = "SELECT minqs,maxqs FROM meta LIMIT 1;";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);


	if (ret != SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		return NULL;
	}

	PYFASTX_SQLITE_CALL(
		minqs = sqlite3_column_int(stmt, 0);
		maxqs = sqlite3_column_int(stmt, 1);
		sqlite3_finalize(stmt);
	);

	platforms = PyList_New(0);

	//check fastq platform
	if (minqs < 33 || maxqs > 126) {
		//return PyErr_Format(PyExc_TypeError, "Quality values corrupt. found [%d, %d] where [33, 104] was expected", minqs, maxqs);
		platform = Py_BuildValue("s", "Unknown");
		PyList_Append(platforms, platform);
		Py_DECREF(platform);
		return platforms;
	}

	if (minqs >= 33 && maxqs <= 73) {
		platform = Py_BuildValue("s", "Sanger Phred+33");
		PyList_Append(platforms, platform);
		Py_DECREF(platform);
	}

	if (minqs >= 33 && maxqs <= 74) {
		platform = Py_BuildValue("s", "Illumina 1.8+ Phred+33");
		PyList_Append(platforms, platform);
		Py_DECREF(platform);
	}

	if (minqs >= 59 && maxqs <= 104) {
		platform = Py_BuildValue("s", "Solexa Solexa+64");
		PyList_Append(platforms, platform);
		Py_DECREF(platform);
	}

	if (minqs >= 64 && maxqs <= 104) {
		platform = Py_BuildValue("s", "Illumina 1.3+ Phred+64");
		PyList_Append(platforms, platform);
		Py_DECREF(platform);
	}

	if (minqs >= 66 && maxqs <= 104) {
		platform = Py_BuildValue("s", "Illumina 1.5+ Phred+64");
		PyList_Append(platforms, platform);
		Py_DECREF(platform);
	}

	if (minqs >= 33 && maxqs <= 126) {
		platform = Py_BuildValue("s", "PacBio HiFi Phred+33");
		PyList_Append(platforms, platform);
		Py_DECREF(platform);
	}

	return platforms;
}

PyObject* pyfastx_fastq_phred(pyfastx_Fastq *self, void* closure) {
	if (!self->middle->phred) {
		pyfastx_fastq_calc_composition(self);
	}

	return Py_BuildValue("i", self->middle->phred);
}

PyObject* pyfastx_fastq_minqual(pyfastx_Fastq *self, void* closure) {
	if (!self->minqual) {
		pyfastx_fastq_calc_composition(self);
	}

	return Py_BuildValue("i", self->minqual);
}

PyObject* pyfastx_fastq_maxqual(pyfastx_Fastq *self, void* closure) {
	if (!self->maxqual) {
		pyfastx_fastq_calc_composition(self);
	}

	return Py_BuildValue("i", self->maxqual);
}

PyObject* pyfastx_fastq_minlen(pyfastx_Fastq *self, void* closure) {
	int ret;
	sqlite3_stmt *stmt;

	if (!self->minlen) {
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index_db, "SELECT minlen FROM meta", -1, &stmt, NULL);
			ret = sqlite3_step(stmt);
		);

		if (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(self->minlen=sqlite3_column_int64(stmt, 0));
		}
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		stmt = NULL;

		if (!self->minlen) {
			PYFASTX_SQLITE_CALL(
				sqlite3_prepare_v2(self->index_db, "SELECT MIN(rlen) FROM read", -1, &stmt, NULL);
				ret = sqlite3_step(stmt);
			);

			if (ret == SQLITE_ROW) {
				PYFASTX_SQLITE_CALL(self->minlen=sqlite3_column_int64(stmt, 0));
			}

			sqlite3_finalize(stmt);
		}

	}

	return Py_BuildValue("n", self->minlen);
}

PyObject* pyfastx_fastq_maxlen(pyfastx_Fastq *self, void* closure) {
	int ret;
	sqlite3_stmt *stmt;

	if (!self->maxlen) {
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index_db, "SELECT maxlen FROM meta", -1, &stmt, NULL);
			ret = sqlite3_step(stmt);
		);

		if (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(self->maxlen=sqlite3_column_int64(stmt, 0));
		}
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		stmt = NULL;

		if (!self->maxlen) {
			PYFASTX_SQLITE_CALL(
				sqlite3_prepare_v2(self->index_db, "SELECT MAX(rlen) FROM read", -1, &stmt, NULL);
				ret = sqlite3_step(stmt);
			);

			if (ret == SQLITE_ROW) {
				PYFASTX_SQLITE_CALL(self->maxlen=sqlite3_column_int64(stmt, 0));
			}

			sqlite3_finalize(stmt);
		}

	}

	return Py_BuildValue("n", self->maxlen);
}

PyObject* pyfastx_fastq_gc_content(pyfastx_Fastq *self, void* closure) {
	int ret;
	const char *sql;

	sqlite3_stmt *stmt;

	Py_ssize_t a, c, g, t;

	if (self->gc_content) {
		return Py_BuildValue("f", self->gc_content);
	}

	pyfastx_fastq_calc_composition(self);

	sql = "SELECT * FROM base LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);
	
	if (ret != SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_SetString(PyExc_RuntimeError, "could not calculate gc content");
		return NULL;
	}

	PYFASTX_SQLITE_CALL(
		a = sqlite3_column_int64(stmt, 0);
		c = sqlite3_column_int64(stmt, 1);
		g = sqlite3_column_int64(stmt, 2);
		t = sqlite3_column_int64(stmt, 3);
		sqlite3_finalize(stmt);
	);

	self->gc_content = (float)(g+c)/(a+c+g+t)*100;

	return Py_BuildValue("f", self->gc_content);
}

PyObject* pyfastx_fastq_composition(pyfastx_Fastq *self, void* closure) {
	int ret;
	const char *sql;

	sqlite3_stmt *stmt;

	Py_ssize_t a, c, g, t, n;

	pyfastx_fastq_calc_composition(self);
	
	sql = "SELECT * FROM base LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);
	
	if (ret != SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_SetString(PyExc_RuntimeError, "could not get composition");
		return NULL;
	}

	PYFASTX_SQLITE_CALL(
		a = sqlite3_column_int64(stmt, 0);
		c = sqlite3_column_int64(stmt, 1);
		g = sqlite3_column_int64(stmt, 2);
		t = sqlite3_column_int64(stmt, 3);
		n = sqlite3_column_int64(stmt, 4);
		sqlite3_finalize(stmt);
	);

	return Py_BuildValue("{s:n,s:n,s:n,s:n,s:n}","A",a,"C",c,"G",g,"T",t,"N",n);
}

PyObject* pyfastx_fastq_isgzip(pyfastx_Fastq *self, void* closure) {
	if (self->middle->gzip_format) {
		Py_RETURN_TRUE;
	}

	Py_RETURN_FALSE;
}

PyObject *pyfastx_fastq_keys(pyfastx_Fastq *self, void* closure) {
	return pyfastx_fastq_keys_create(self->index_db, self->read_counts);
}

static PySequenceMethods pyfastx_fastq_as_sequence = {
	.sq_contains = (objobjproc)pyfastx_fastq_contains,
};

static PyMappingMethods pyfastx_fastq_as_mapping = {
	.mp_length = (lenfunc)pyfastx_fastq_length,
	.mp_subscript = (binaryfunc)pyfastx_fastq_subscript,
};

static PyMethodDef pyfastx_fastq_methods[] = {
	{"build_index", (PyCFunction)pyfastx_fastq_build_index, METH_NOARGS, NULL},
	{"keys", (PyCFunction)pyfastx_fastq_keys, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyGetSetDef pyfastx_fastq_getsets[] = {
	{"encoding_type", (getter)pyfastx_fastq_guess_encoding_type, NULL, NULL, NULL},
	{"phred", (getter)pyfastx_fastq_phred, NULL, NULL, NULL},
	{"gc_content", (getter)pyfastx_fastq_gc_content, NULL, NULL, NULL},
	{"composition", (getter)pyfastx_fastq_composition, NULL, NULL, NULL},
	{"maxlen", (getter)pyfastx_fastq_maxlen, NULL, NULL, NULL},
	{"minlen", (getter)pyfastx_fastq_minlen, NULL, NULL, NULL},
	{"maxqual", (getter)pyfastx_fastq_maxqual, NULL, NULL, NULL},
	{"minqual", (getter)pyfastx_fastq_minqual, NULL, NULL, NULL},
	{"is_gzip", (getter)pyfastx_fastq_isgzip, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef pyfastx_fastq_members[] = {
	{"file_name", T_OBJECT, offsetof(pyfastx_Fastq, file_obj), READONLY},
	{"size", T_PYSSIZET, offsetof(pyfastx_Fastq, seq_length), READONLY},
	{"avglen", T_DOUBLE, offsetof(pyfastx_Fastq, avg_length), READONLY},
	{NULL}
};

PyTypeObject pyfastx_FastqType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Fastq",
    .tp_basicsize = sizeof(pyfastx_Fastq),
    .tp_dealloc = (destructor)pyfastx_fastq_dealloc,
    .tp_repr = (reprfunc)pyfastx_fastq_repr,
    .tp_as_sequence = &pyfastx_fastq_as_sequence,
    .tp_as_mapping = &pyfastx_fastq_as_mapping,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_iter = (getiterfunc)pyfastx_fastq_iter,
    .tp_iternext = (iternextfunc)pyfastx_fastq_next,
    .tp_methods = pyfastx_fastq_methods,
    .tp_members = pyfastx_fastq_members,
    .tp_getset = pyfastx_fastq_getsets,
    .tp_new = pyfastx_fastq_new,
};