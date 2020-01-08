#include "fastq.h"
#include "read.h"
#include "structmember.h"

void pyfastx_fastq_build_index(pyfastx_Fastq *self) {
	sqlite3_stmt *stmt;
	int64_t soff = 0, qoff = 0, pos = 0;
	uint64_t size = 0;
	char* name = NULL;
	char* space;
	int l, j, rlen = 0;
	int dlen = 0;
	const char *sql;
	kstring_t line = {0, 0, 0};
	uint64_t line_num = 0;

	sql = " \
		CREATE TABLE read ( \
			ID INTEGER PRIMARY KEY, --read id \n \
			name TEXT, --read name \n \
			dlen INTEGER, --description length \n \
			rlen INTEGER, --read length \n \
			soff INTEGER, --read seq offset \n \
			qoff INTEGER --read qual offset \n \
		); \
		CREATE TABLE meta ( \
			count INTEGER, --read count \n \
			size INTEGER --all read length \n \
		); \
		CREATE TABLE gzindex (  \
			ID INTEGER PRIMARY KEY,  \
			content BLOB  \
		); \
		CREATE TABLE base ( \
			a INTEGER,  \
			c INTEGER,  \
			g INTEGER,  \
			t INTEGER,  \
			n INTEGER  \
		); \
		CREATE TABLE qual ( \
			minqs INTEGER, --max quality score \n \
			maxqs INTEGER, --min quality score \n \
			phred INTEGER --phred value \n \
		);";


	if(sqlite3_open(self->index_file, &self->index_db) != SQLITE_OK){
		PyErr_SetString(PyExc_ConnectionError, sqlite3_errmsg(self->index_db));
		return;
	}

	if (sqlite3_exec(self->index_db, sql, NULL, NULL, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return;
	}

	sql = "PRAGMA synchronous = OFF; BEGIN TRANSACTION;";

	if(sqlite3_exec(self->index_db, sql, NULL, NULL, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return;
	}

	sql = "INSERT INTO read VALUES (?,?,?,?,?,?);";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	
	Py_BEGIN_ALLOW_THREADS

	gzrewind(self->gzfd);
	ks_rewind(self->ks);

	while ((l=ks_getuntil(self->ks, '\n', &line, 0)) >= 0) {
		++line_num;
		++l;

		j = line_num % 4;

		switch(j) {
			case 1:
				name = (char *)malloc(line.l);
				memcpy(name, line.s+1, line.l);

				dlen = line.l;

				if ((space = strchr(name, ' ')) != NULL) {
					*space = '\0';
				}

				if ((space = strchr(name, '\r')) != NULL) {
					*space = '\0';
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
				sqlite3_bind_text(stmt, 2, name, -1, NULL);
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
	
	sqlite3_exec(self->index_db, "CREATE INDEX readidx ON read (name);", NULL, NULL, NULL);
	sqlite3_exec(self->index_db, "COMMIT;", NULL, NULL, NULL);

	self->read_counts = line_num/4;
	self->seq_length = size;
	sql = "INSERT INTO meta VALUES (?,?);";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int64(stmt, 1, self->read_counts);
	sqlite3_bind_int64(stmt, 2, self->seq_length);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	Py_END_ALLOW_THREADS

	//if is gzip build gzip index
	if (self->gzip_format) {
		pyfastx_build_gzip_index(self->gzip_index, self->index_db, self->index_file);
	}

	free(line.s);
}

void pyfastx_fastq_load_index(pyfastx_Fastq *self) {
	sqlite3_stmt* stmt;
	const char* sql;
	
	if(sqlite3_open(self->index_file, &self->index_db) != SQLITE_OK){
		PyErr_SetString(PyExc_ConnectionError, sqlite3_errmsg(self->index_db));
		return;
	}

	//calculate attributes
	sql = "SELECT * FROM meta LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	if (sqlite3_step(stmt) == SQLITE_ROW) {
		self->read_counts = sqlite3_column_int64(stmt, 0);
		self->seq_length = sqlite3_column_int64(stmt, 1);
	} else {
		PyErr_SetString(PyExc_RuntimeError, "Can not get read counts and length");
		return;
	}
	sqlite3_reset(stmt);

	sql = "SELECT phred FROM qual LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	if (sqlite3_step(stmt) == SQLITE_ROW) {
		self->phred = sqlite3_column_int(stmt, 0);
	} else {
		self->phred = 0;
	}

	if(self->gzip_format){
		pyfastx_load_gzip_index(self->gzip_index, self->index_db, self->index_file);
	}
}

PyObject *pyfastx_fastq_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	char *file_name;
	int phred = 0;
	int build_index = 1;
	int composition = 0;

	static char* keywords[] = {"file_name", "phred", "build_index", "composition", NULL};

	pyfastx_Fastq *obj;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|iii", keywords, &file_name, &phred, &build_index, &composition)) {
		return NULL;
	}

	//check input sequence file is whether exists
	if (!file_exists(file_name)) {
		PyErr_Format(PyExc_FileExistsError, "input fastq file %s does not exists", file_name);
		return NULL;
	}

	obj = (pyfastx_Fastq *)type->tp_alloc(type, 0);
	if (!obj) {
		return NULL;
	}

	obj->file_name = (char *)malloc(strlen(file_name) + 1);
	strcpy(obj->file_name, file_name);

	//check input file is gzip or not
	obj->gzip_format = is_gzip_format(file_name);

	//initial kstream and kseq
	obj->gzfd = gzopen(file_name, "rb");
	obj->ks = ks_init(obj->gzfd);
	obj->kseq = kseq_init(obj->gzfd);

	//create index file
	obj->index_file = (char *)malloc(strlen(file_name) + 5);
	strcpy(obj->index_file, file_name);
	strcat(obj->index_file, ".fxi");

	obj->fd = fopen(file_name, "rb");

	//initail index connection
	obj->index_db = NULL;

	obj->has_index = build_index;

	//inital phred
	obj->phred = phred;
	obj->gc_content = 0;

	//is gzip file
	if (obj->gzip_format) {
		obj->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
		zran_init(obj->gzip_index, obj->fd, 4194304, 32768, 1048576, ZRAN_AUTO_BUILD);
	}

	if (file_exists(obj->index_file)) {
		pyfastx_fastq_load_index(obj);
	} else if (build_index) {
		pyfastx_fastq_build_index(obj);
	}

	if (build_index && composition) {
		pyfastx_fastq_calc_composition(obj);
	}

	return (PyObject *)obj;
}

void pyfastx_fastq_dealloc(pyfastx_Fastq *self) {
	if (self->index_db) {
		sqlite3_close(self->index_db);
	}

	if (self->gzip_format) {
		zran_free(self->gzip_index);
	}

	ks_destroy(self->ks);
	kseq_destroy(self->kseq);
	fclose(self->fd);
	gzclose(self->gzfd);

	Py_TYPE(self)->tp_free((PyObject *)self);
}

int pyfastx_fastq_length(pyfastx_Fastq *self) {
	return self->read_counts;
}

PyObject* pyfastx_fastq_make_read(pyfastx_Fastq *self, sqlite3_stmt *stmt) {
	//int a, c, g, t, n;
	int nbytes;

	pyfastx_Read *read = PyObject_New(pyfastx_Read, &pyfastx_ReadType);

	if (!read) {
		return NULL;
	}

	read->id = sqlite3_column_int64(stmt, 0);
	nbytes = sqlite3_column_bytes(stmt, 1);
	read->name = (char *)malloc(nbytes + 1);
	memcpy(read->name, (char *)sqlite3_column_text(stmt, 1), nbytes);
	read->name[nbytes] = '\0';
	read->desc_len = sqlite3_column_int(stmt, 2);
	read->read_len = sqlite3_column_int(stmt, 3);
	read->seq_offset = sqlite3_column_int64(stmt, 4);
	read->qual_offset = sqlite3_column_int64(stmt, 5);
	read->gzfd = self->gzfd;
	read->gzip_index = self->gzip_index;
	read->gzip_format = self->gzip_format;
	read->phred = self->phred;
	read->seq = NULL;
	read->qual = NULL;

	sqlite3_finalize(stmt);

	Py_INCREF(read);
	return (PyObject *)read;
}

PyObject* pyfastx_fastq_get_read_by_id(pyfastx_Fastq *self, uint32_t read_id) {
	sqlite3_stmt *stmt;

	const char* sql = "SELECT * FROM read WHERE ID=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, read_id);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		PyErr_SetString(PyExc_IndexError, "Index Error");
		return NULL;
	}

	return pyfastx_fastq_make_read(self, stmt);
}

PyObject* pyfastx_fastq_get_read_by_name(pyfastx_Fastq *self, char* name) {
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	
	//select sql statement, chrom indicates seq name or chromomsome
	const char* sql = "SELECT * FROM read WHERE name=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, name, -1, NULL);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		PyErr_SetString(PyExc_KeyError, name);
		return NULL;
	}

	return pyfastx_fastq_make_read(self, stmt);
}

PyObject* pyfastx_fastq_subscript(pyfastx_Fastq *self, PyObject *item) {
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
		char *key = PyUnicode_AsUTF8(item);

		return pyfastx_fastq_get_read_by_name(self, key);

	} else {
		PyErr_SetString(PyExc_KeyError, "key error");
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

	if (!PyUnicode_Check(key)) {
		return 0;
	}

	name = PyUnicode_AsUTF8(key);
	sql = "SELECT * FROM read WHERE name=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, name, -1, NULL);

	if (sqlite3_step(stmt) == SQLITE_ROW) {
		return 1;
	}

	return 0;
}

PyObject *pyfastx_fastq_iter(pyfastx_Fastq *self) {
	gzrewind(self->gzfd);
	kseq_rewind(self->kseq);

	if (self->has_index) {
		self->iter_id = 0;
	}
	
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fastq_next(pyfastx_Fastq *self) {
	if (self->has_index) {
		++self->iter_id;

		if (self->iter_id <= self->read_counts) {
			return pyfastx_fastq_get_read_by_id(self, self->iter_id);
		}
	} else {
		if (kseq_read(self->kseq) >= 0) {
			return Py_BuildValue("sss", self->kseq->name.s, self->kseq->seq.s, self->kseq->qual.s);
		}
	}

	return NULL;
}

void pyfastx_fastq_calc_composition(pyfastx_Fastq *self) {
	sqlite3_stmt *stmt;
	kstring_t line = {0, 0, 0};
	uint64_t a = 0, c = 0, g = 0, t = 0, n = 0;
	uint64_t line_num = 0;
	int i, j;

	//min and max quality score
	int minqs = 104, maxqs = 33, phred = 0;

	const char *sql = "SELECT * FROM base LIMIT 1";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	if (sqlite3_step(stmt) == SQLITE_ROW) {
		return;
	}
	
	Py_BEGIN_ALLOW_THREADS

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
					continue;
				}

				if (line.s[i] < minqs) {
					minqs = line.s[i];
				}

				if (line.s[i] > maxqs) {
					maxqs = line.s[i];
				}
			}
		}
	}

	sql = "INSERT INTO base VALUES (?,?,?,?,?);";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int64(stmt, 1, a);
	sqlite3_bind_int64(stmt, 2, c);
	sqlite3_bind_int64(stmt, 3, g);
	sqlite3_bind_int64(stmt, 4, t);
	sqlite3_bind_int64(stmt, 5, n);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	if (phred == 0) {
		if (maxqs > 74) {
			phred = 64;
		}

		if (minqs < 59) {
			phred = 33;
		}
	}

	//insert platform into index file
	sql = "INSERT INTO qual VALUES (?,?,?);";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, minqs);
	sqlite3_bind_int(stmt, 2, maxqs);
	sqlite3_bind_int(stmt, 3, phred);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	Py_END_ALLOW_THREADS
}

PyObject* pyfastx_fastq_guess_encoding_type(pyfastx_Fastq* self, void* closure) {
	sqlite3_stmt *stmt;
	int maxqs, minqs;
	const char *sql;
	PyObject* platforms;
	PyObject* platform;

	pyfastx_fastq_calc_composition(self);

	sql = "SELECT * FROM qual LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	if (sqlite3_step(stmt) != SQLITE_ROW) {
		return NULL;
	}

	minqs = sqlite3_column_int(stmt, 0);
	maxqs = sqlite3_column_int(stmt, 1);

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
	if (self->phred == 0) {
		sqlite3_stmt *stmt;
		const char *sql;
		
		pyfastx_fastq_calc_composition(self);

		sql = "SELECT phred FROM qual LIMIT 1;";
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		if (sqlite3_step(stmt) != SQLITE_ROW) {
			return NULL;
		}
		self->phred = sqlite3_column_int(stmt, 0);
	}

	return Py_BuildValue("i", self->phred);
}

PyObject* pyfastx_fastq_gc_content(pyfastx_Fastq *self, void* closure) {
	sqlite3_stmt *stmt;
	int64_t a, c, g, t;
	const char *sql;

	if (self->gc_content > 0) {
		return Py_BuildValue("f", self->gc_content);
	}

	pyfastx_fastq_calc_composition(self);

	sql = "SELECT * FROM base LIMIT 1";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	
	if (sqlite3_step(stmt) != SQLITE_ROW) {
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return NULL;
	}

	a = sqlite3_column_int64(stmt, 0);
	c = sqlite3_column_int64(stmt, 1);
	g = sqlite3_column_int64(stmt, 2);
	t = sqlite3_column_int64(stmt, 3);

	self->gc_content = (float)(g+c)/(a+c+g+t)*100;

	return Py_BuildValue("f", self->gc_content);
}

PyObject* pyfastx_fastq_composition(pyfastx_Fastq *self, void* closure) {
	sqlite3_stmt *stmt;
	int64_t a, c, g, t, n;
	const char *sql;

	pyfastx_fastq_calc_composition(self);
	
	sql = "SELECT * FROM base LIMIT 1";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	
	if (sqlite3_step(stmt) != SQLITE_ROW) {
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return NULL;
	}

	a = sqlite3_column_int64(stmt, 0);
	c = sqlite3_column_int64(stmt, 1);
	g = sqlite3_column_int64(stmt, 2);
	t = sqlite3_column_int64(stmt, 3);
	n = sqlite3_column_int64(stmt, 4);

	return Py_BuildValue("{s:I,s:I,s:I,s:I,s:I}","A",a,"C",c,"G",g,"T",t,"N",n);
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

static PyGetSetDef pyfastx_fastq_getsets[] = {
	{"encoding_type", (getter)pyfastx_fastq_guess_encoding_type, NULL, NULL, NULL},
	{"phred", (getter)pyfastx_fastq_phred, NULL, NULL, NULL},
	{"gc_content", (getter)pyfastx_fastq_gc_content, NULL, NULL, NULL},
	{"composition", (getter)pyfastx_fastq_composition, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef pyfastx_fastq_members[] = {
	{"file_name", T_STRING, offsetof(pyfastx_Fastq, file_name), READONLY},
	{"size", T_LONG, offsetof(pyfastx_Fastq, seq_length), READONLY},
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
    0,          /* tp_methods */
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