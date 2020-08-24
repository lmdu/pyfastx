#include "fastq.h"
#include "read.h"
#include "structmember.h"

void pyfastx_fastq_create_index(pyfastx_Fastq *self) {
	sqlite3_stmt *stmt;
	int64_t soff = 0, qoff = 0, pos = 0;
	uint64_t size = 0;
	kstring_t name = {0, 0, 0};
	char* space;
	int l, j, rlen = 0;
	int dlen = 0;
	const char *sql;
	kstring_t line = {0, 0, 0};
	uint64_t line_num = 0;
	int ret;

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
		PyErr_Format(PyExc_ConnectionError, "can not open index file %s", self->index_file);
		return;
	}

	PYFASTX_SQLITE_CALL(ret=sqlite3_exec(self->index_db, sql, NULL, NULL, NULL));
	if (ret != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, "can not create index table");
		return;
	}

	sql = "PRAGMA synchronous = OFF; BEGIN TRANSACTION;";
	PYFASTX_SQLITE_CALL(ret=sqlite3_exec(self->index_db, sql, NULL, NULL, NULL));
	if(ret != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, "can not begin transaction");
		return;
	}
	
	sql = "INSERT INTO read VALUES (?,?,?,?,?,?);";
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));

	gzrewind(self->gzfd);
	ks_rewind(self->ks);

	Py_BEGIN_ALLOW_THREADS

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
				sqlite3_bind_null(stmt, 1);
				sqlite3_bind_text(stmt, 2, name.s, name.l, SQLITE_STATIC);
				sqlite3_bind_int(stmt, 3, dlen);
				sqlite3_bind_int(stmt, 4, rlen);
				sqlite3_bind_int64(stmt, 5, soff);
				sqlite3_bind_int64(stmt, 6, qoff);
				sqlite3_step(stmt);
				sqlite3_reset(stmt);
				break;
		}
		pos += l;
	}

	sqlite3_finalize(stmt);
	stmt = NULL;

	sqlite3_exec(self->index_db, "CREATE INDEX readidx ON read (name);", NULL, NULL, NULL);
	sqlite3_exec(self->index_db, "COMMIT;", NULL, NULL, NULL);

	self->read_counts = line_num/4;
	self->seq_length = size;
	self->avg_length = size*1.0/self->read_counts;
	sql = "INSERT INTO stat VALUES (?,?,?);";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int64(stmt, 1, self->read_counts);
	sqlite3_bind_int64(stmt, 2, self->seq_length);
	sqlite3_bind_double(stmt, 3, self->avg_length);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	Py_END_ALLOW_THREADS

	free(line.s);
	free(name.s);

	//if is gzip build gzip index
	if (self->gzip_format) {
		pyfastx_build_gzip_index(self->gzip_index, self->index_db);
	}
}

void pyfastx_fastq_load_index(pyfastx_Fastq *self) {
	sqlite3_stmt* stmt;
	const char* sql;
	int ret;

	PYFASTX_SQLITE_CALL(ret=sqlite3_open(self->index_file, &self->index_db));

	if (ret != SQLITE_OK){
		PyErr_Format(PyExc_ConnectionError, "can not open index file %s", self->index_file);
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
		PYFASTX_SQLITE_CALL(self->phred = sqlite3_column_int(stmt, 0));
	} else {
		self->phred = 0;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	if(self->gzip_format){
		pyfastx_load_gzip_index(self->gzip_index, self->index_db);
	}
}

PyObject *pyfastx_fastq_build_index(pyfastx_Fastq *self){
	if (file_exists(self->index_file)) {
		pyfastx_fastq_load_index(self);
	} else {
		pyfastx_fastq_create_index(self);
	}

	Py_RETURN_TRUE;
}

PyObject *pyfastx_fastq_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	PyObject *file_obj;
	char *file_name;
	Py_ssize_t file_len;

	int phred = 0;
	int build_index = 1;
	int full_index = 0;

	static char* keywords[] = {"file_name", "phred", "build_index", "full_index", NULL};

	pyfastx_Fastq *obj;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|iii", keywords, &file_obj, &phred, &build_index, &full_index)) {
		return NULL;
	}

	//check input sequence file is whether exists
	file_name = (char *)PyUnicode_AsUTF8AndSize(file_obj, &file_len);
	
	if (!file_name) {
		PyErr_Format(PyExc_ValueError, "the input file name is not a right string");
		return NULL;
	}

	if (!file_exists(file_name)) {
		PyErr_Format(PyExc_FileExistsError, "input fastq file %s does not exists", file_name);
		return NULL;
	}

	obj = (pyfastx_Fastq *)type->tp_alloc(type, 0);
	if (!obj) {
		return NULL;
	}

	obj->file_name = (char *)malloc((int)file_len + 1);
	strcpy(obj->file_name, file_name);

	//check input file is gzip or not
	obj->gzip_format = is_gzip_format(file_name);

	//initial kstream and kseq
	obj->gzfd = gzopen(file_name, "rb");
	obj->ks = ks_init(obj->gzfd);
	obj->kseq = kseq_init(obj->gzfd);

	//check is correct fastq format
	if (!fastq_validator(obj->gzfd)) {
		PyErr_Format(PyExc_RuntimeError, "%s is not plain or gzip compressed fastq formatted file", file_name);
		return NULL;
	}

	//create index file
	obj->index_file = (char *)malloc((int)file_len + 5);
	strcpy(obj->index_file, file_name);
	strcat(obj->index_file, ".fxi");

	obj->fd = fopen(file_name, "rb");

	//initail index connection
	obj->index_db = 0;
	obj->iter_stmt = NULL;

	obj->has_index = build_index;

	//initialize attribute
	obj->phred = phred;
	obj->gc_content = 0;
	obj->minlen = 0;
	obj->maxlen = 0;
	obj->minqual = 0;
	obj->maxqual = 0;

	//is gzip file
	if (obj->gzip_format) {
		obj->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
		zran_init(obj->gzip_index, obj->fd, 1048576, 32768, 16384, ZRAN_AUTO_BUILD);
	}

	if (file_exists(obj->index_file)) {
		pyfastx_fastq_load_index(obj);
	} else if (build_index) {
		pyfastx_fastq_create_index(obj);
	}

	if (build_index && full_index) {
		pyfastx_fastq_calc_composition(obj);
	}

	//initialize cache buffer
	obj->cache_buff = NULL;
	obj->cache_soff = 0;
	obj->cache_eoff = 0;

	return (PyObject *)obj;
}

void pyfastx_fastq_dealloc(pyfastx_Fastq *self) {
	if (self->iter_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->iter_stmt));
	}

	if (self->index_db) {
		PYFASTX_SQLITE_CALL(sqlite3_close(self->index_db));
	}

	if (self->gzip_format) {
		zran_free(self->gzip_index);
	}

	if (self->cache_buff) {
		free(self->cache_buff);
	}

	ks_destroy(self->ks);
	kseq_destroy(self->kseq);
	fclose(self->fd);
	gzclose(self->gzfd);

	Py_TYPE(self)->tp_free((PyObject *)self);
}

uint64_t pyfastx_fastq_length(pyfastx_Fastq *self) {
	return self->read_counts;
}

PyObject* pyfastx_fastq_make_read(pyfastx_Fastq *self, sqlite3_stmt *stmt) {
	int nbytes;

	//pyfastx_Read *read = PyObject_New(pyfastx_Read, &pyfastx_ReadType);
	pyfastx_Read *read = (pyfastx_Read *)PyObject_CallObject((PyObject *)&pyfastx_ReadType, NULL);

	PYFASTX_SQLITE_CALL(
		read->id = sqlite3_column_int64(stmt, 0);
		nbytes = sqlite3_column_bytes(stmt, 1);
		read->name = (char *)malloc(nbytes + 1);
		memcpy(read->name, sqlite3_column_text(stmt, 1), nbytes);
		read->name[nbytes] = '\0';
		read->desc_len = sqlite3_column_int(stmt, 2);
		read->read_len = sqlite3_column_int(stmt, 3);
		read->seq_offset = sqlite3_column_int64(stmt, 4);
		read->qual_offset = sqlite3_column_int64(stmt, 5);
		//sqlite3_finalize(stmt);
	);

	read->fastq = self;
	read->seq = NULL;
	read->qual = NULL;
	read->raw = NULL;
	read->desc = NULL;

	//Py_INCREF(read);
	return (PyObject *)read;
}

PyObject* pyfastx_fastq_get_read_by_id(pyfastx_Fastq *self, uint64_t read_id) {
	sqlite3_stmt *stmt;
	PyObject *obj;
	int ret;

	const char* sql = "SELECT * FROM read WHERE ID=? LIMIT 1;";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int(stmt, 1, read_id);
		ret = sqlite3_step(stmt);
	);

	if (ret != SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_SetString(PyExc_IndexError, "Index Error");
		return NULL;
	}

	obj = pyfastx_fastq_make_read(self, stmt);
	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return obj;
}

PyObject* pyfastx_fastq_get_read_by_name(pyfastx_Fastq *self, char* name) {
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	PyObject *obj;
	int ret;
	
	//select sql statement, chrom indicates seq name or chromomsome
	const char* sql = "SELECT * FROM read WHERE name=? LIMIT 1;";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret != SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_Format(PyExc_KeyError, "%s does not exist in fastq file", name);
		return NULL;
	}

	obj = pyfastx_fastq_make_read(self, stmt);
	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return obj;
}

PyObject* pyfastx_fastq_subscript(pyfastx_Fastq *self, PyObject *item) {
	self->iterating = 0;

	if (PyIndex_Check(item)) {
		Py_ssize_t i;
		i = PyNumber_AsSsize_t(item, PyExc_IndexError);

		if (i < 0) {
			i += self->read_counts;
		}

		if(i >= self->read_counts){
			PyErr_SetString(PyExc_IndexError, "index out of range");
			return NULL;
		}

		return pyfastx_fastq_get_read_by_id(self, i+1);
		
	} else if (PyUnicode_Check(item)) {
		return pyfastx_fastq_get_read_by_name(self, (char *)PyUnicode_AsUTF8(item));

	} else {
		PyErr_SetString(PyExc_KeyError, "the key must be index number or read name");
		return NULL;
	}
}

PyObject* pyfastx_fastq_repr(pyfastx_Fastq *self) {
	return PyUnicode_FromFormat("<Fastq> %s contains %ld reads", self->file_name, self->read_counts);
}

int pyfastx_fastq_contains(pyfastx_Fastq *self, PyObject *key) {
	sqlite3_stmt *stmt;
	char *name;
	const char *sql;
	int ret;

	if (!PyUnicode_Check(key)) {
		return 0;
	}

	name = (char *)PyUnicode_AsUTF8(key);
	sql = "SELECT * FROM read WHERE name=? LIMIT 1;";

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	return ret==SQLITE_ROW ? 1 : 0;
}

PyObject *pyfastx_fastq_iter(pyfastx_Fastq *self) {
	gzrewind(self->gzfd);
	rewind(self->fd);
	
	if (self->has_index) {
		self->iterating = 1;

		if (!self->cache_buff) {
			self->cache_buff = (char *)malloc(1048576);
		}

		PYFASTX_SQLITE_CALL(
			sqlite3_finalize(self->iter_stmt);
			self->iter_stmt = NULL;
			sqlite3_prepare_v2(self->index_db, "SELECT * FROM read", -1, &self->iter_stmt, NULL);
		);
	} else {
		kseq_rewind(self->kseq);
	}
	
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fastq_next(pyfastx_Fastq *self) {
	if (self->has_index) {
		int ret;
		PYFASTX_SQLITE_CALL(ret = sqlite3_step(self->iter_stmt));

		if (ret == SQLITE_ROW) {
			return pyfastx_fastq_make_read(self, self->iter_stmt);
		}

	} else {
		if (kseq_read(self->kseq) >= 0) {
			return Py_BuildValue("sss", self->kseq->name.s, self->kseq->seq.s, self->kseq->qual.s);
		}
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(self->iter_stmt));
	self->iter_stmt = NULL;
	self->iterating = 0;

	if (self->cache_buff) {
		free(self->cache_buff);
		self->cache_buff = NULL;
	}

	return NULL;
}

void pyfastx_fastq_calc_composition(pyfastx_Fastq *self) {
	sqlite3_stmt *stmt;
	kstring_t line = {0, 0, 0};
	uint64_t a = 0, c = 0, g = 0, t = 0, n = 0;
	uint64_t line_num = 0;
	int i, j;
	int ret;

	//min and max read length
	int minlen=10000000, maxlen=0;

	//min and max quality score
	int minqs = 104, maxqs = 33, phred = 0;

	const char *sql = "SELECT * FROM meta LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			if (!self->maxlen)
				self->maxlen = sqlite3_column_int(stmt, 0);
			if (!self->minlen)
				self->minlen = sqlite3_column_int(stmt, 1);
			if (!self->minqual)
				self->minqual = sqlite3_column_int(stmt, 2);
			if (!self->maxqual)
				self->maxqual = sqlite3_column_int(stmt, 3);
			if (!self->phred)
				self->phred = sqlite3_column_int(stmt, 4);
		);

		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		return;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
	stmt = NULL;

	gzrewind(self->gzfd);
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
		sqlite3_bind_int(stmt, 1, maxlen);
		sqlite3_bind_int(stmt, 2, minlen);
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
	self->phred = phred;
}

PyObject* pyfastx_fastq_guess_encoding_type(pyfastx_Fastq* self, void* closure) {
	sqlite3_stmt *stmt;
	int maxqs, minqs;
	const char *sql;
	PyObject* platforms;
	PyObject* platform;
	int ret;

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
	if (minqs < 33 || maxqs > 104) {
		return PyErr_Format(PyExc_TypeError, "Quality values corrupt. found [%d, %d] where [33, 104] was expected", minqs, maxqs);
	}

	if (minqs >= 33 && maxqs <= 73) {
		platform = Py_BuildValue("s", "Sanger Phred+33");
		PyList_Append(platforms, platform);
	}

	if (minqs >= 33 && maxqs <= 74) {
		platform = Py_BuildValue("s", "Illumina 1.8+ Phred+33");
		PyList_Append(platforms, platform);
	}

	if (minqs >= 59 && maxqs <= 104) {
		platform = Py_BuildValue("s", "Solexa Solexa+64");
		PyList_Append(platforms, platform);
	}

	if (minqs >= 64 && maxqs <= 104) {
		platform = Py_BuildValue("s", "Illumina 1.3+ Phred+64");
		PyList_Append(platforms, platform);
	}

	if (minqs >= 66 && maxqs <= 104) {
		platform = Py_BuildValue("s", "Illumina 1.5+ Phred+64");
		PyList_Append(platforms, platform);
	}

	return platforms;
}

PyObject* pyfastx_fastq_phred(pyfastx_Fastq *self, void* closure) {
	if (!self->phred) {
		pyfastx_fastq_calc_composition(self);
	}

	return Py_BuildValue("i", self->phred);
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
	sqlite3_stmt *stmt;
	int ret;

	if (!self->minlen) {
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index_db, "SELECT minlen FROM meta", -1, &stmt, NULL);
			ret = sqlite3_step(stmt);
		);

		if (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(self->minlen=sqlite3_column_double(stmt, 0));
		}
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		stmt = NULL;

		if (!self->minlen) {
			PYFASTX_SQLITE_CALL(
				sqlite3_prepare_v2(self->index_db, "SELECT min(rlen) FROM read", -1, &stmt, NULL);
				ret = sqlite3_step(stmt);
			);

			if (ret == SQLITE_ROW) {
				PYFASTX_SQLITE_CALL(self->minlen=sqlite3_column_double(stmt, 0));
			}

			sqlite3_finalize(stmt);
		}

	}

	return Py_BuildValue("i", self->minlen);
}

PyObject* pyfastx_fastq_maxlen(pyfastx_Fastq *self, void* closure) {
	sqlite3_stmt *stmt;
	int ret;

	if (!self->maxlen) {
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index_db, "SELECT maxlen FROM meta", -1, &stmt, NULL);
			ret = sqlite3_step(stmt);
		);

		if (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(self->maxlen=sqlite3_column_double(stmt, 0));
		}
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		stmt = NULL;

		if (!self->maxlen) {
			PYFASTX_SQLITE_CALL(
				sqlite3_prepare_v2(self->index_db, "SELECT max(rlen) FROM read", -1, &stmt, NULL);
				ret = sqlite3_step(stmt);
			);

			if (ret == SQLITE_ROW) {
				PYFASTX_SQLITE_CALL(self->maxlen=sqlite3_column_double(stmt, 0));
			}

			sqlite3_finalize(stmt);
		}

	}

	return Py_BuildValue("i", self->maxlen);
}

PyObject* pyfastx_fastq_gc_content(pyfastx_Fastq *self, void* closure) {
	sqlite3_stmt *stmt;
	int64_t a, c, g, t;
	const char *sql;
	int ret;

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
		PyErr_SetString(PyExc_RuntimeError, "can not calculate gc content");
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
	sqlite3_stmt *stmt;
	int64_t a, c, g, t, n;
	const char *sql;
	int ret;

	pyfastx_fastq_calc_composition(self);
	
	sql = "SELECT * FROM base LIMIT 1";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);
	
	if (ret != SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_SetString(PyExc_RuntimeError, "can not get composition");
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

	return Py_BuildValue("{s:L,s:L,s:L,s:L,s:L}","A",a,"C",c,"G",g,"T",t,"N",n);
}

static PySequenceMethods pyfastx_fastq_as_sequence = {
	0, /*sq_length*/
	0, /*sq_concat*/
	0, /*sq_repeat*/
	0, /*sq_item*/
	0, /*sq_slice */
	0, /*sq_ass_item*/
	0, /*sq_ass_splice*/
	(objobjproc)pyfastx_fastq_contains, /*sq_contains*/
	0, /*sq_inplace_concat*/
	0, /*sq_inplace_repeat*/
};

static PyMappingMethods pyfastx_fastq_as_mapping = {
	(lenfunc)pyfastx_fastq_length,
	(binaryfunc)pyfastx_fastq_subscript,
	0,
};

static PyMethodDef pyfastx_fastq_methods[] = {
	{"build_index", (PyCFunction)pyfastx_fastq_build_index, METH_NOARGS, NULL},
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
	{NULL}
};

static PyMemberDef pyfastx_fastq_members[] = {
	{"file_name", T_STRING, offsetof(pyfastx_Fastq, file_name), READONLY},
	{"size", T_ULONGLONG, offsetof(pyfastx_Fastq, seq_length), READONLY},
	{"avglen", T_FLOAT, offsetof(pyfastx_Fastq, avg_length), READONLY},
	{"is_gzip", T_BOOL, offsetof(pyfastx_Fastq, gzip_format), READONLY},
	//{"gc_content", T_FLOAT, offsetof(pyfastx_Fastq, gc_content), READONLY},
	//{"composition", T_OBJECT, offsetof(pyfastx_Fastq, composition), READONLY},
	{NULL}
};

PyTypeObject pyfastx_FastqType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Fastq",                        /* tp_name */
    sizeof(pyfastx_Fastq),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_fastq_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)pyfastx_fastq_repr,                              /* tp_repr */
    0,                              /* tp_as_number */
    &pyfastx_fastq_as_sequence,                   /* tp_as_sequence */
    &pyfastx_fastq_as_mapping,                   /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,             /* tp_flags */
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    (getiterfunc)pyfastx_fastq_iter,     /* tp_iter */
    (iternextfunc)pyfastx_fastq_next,    /* tp_iternext */
    pyfastx_fastq_methods,          /* tp_methods */
    pyfastx_fastq_members,          /* tp_members */
    pyfastx_fastq_getsets,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    pyfastx_fastq_new,              /* tp_new */
};