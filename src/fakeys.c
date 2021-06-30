#include "fakeys.h"
#include "util.h"

char SORTS[][6] = {"ID", "chrom", "slen"};
char ORDERS[][5] = {"ASC", "DESC"};

void pyfastx_fasta_keys_prepare(pyfastx_FastaKeys *self) {
	char *iter_sql;
	char *item_sql;
	char *in_sql;

	if (self->iter_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->iter_stmt));
		self->iter_stmt = NULL;
	}

	if (self->item_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->item_stmt));
		self->item_stmt = NULL;
	}

	if (self->in_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->in_stmt));
		self->in_stmt = NULL;
	}

	iter_sql = sqlite3_mprintf("SELECT chrom FROM seq %s %s %s",
		self->filter ? "WHERE" : "",
		self->filter ? self->filter : "",
		self->order ? self->order : "ORDER BY ID"
	);

	if (self->filter || self->order) {
		item_sql = sqlite3_mprintf("SELECT chrom FROM seq %s %s %s LIMIT 1 OFFSET ?",
			self->filter ? "WHERE" : "",
			self->filter ? self->filter : "",
			self->order ? self->order : "ORDER BY ID"
		);
	} else {
		item_sql = sqlite3_mprintf("SELECT chrom FROM seq WHERE ID=?");
	}

	if (self->filter) {
		in_sql = sqlite3_mprintf("SELECT 1 FROM seq %s %s AND chrom=? LIMIT 1",
			self->filter ? "WHERE" : "",
			self->filter ? self->filter : ""
		);
	} else {
		in_sql = sqlite3_mprintf("SELECT 1 FROM seq WHERE chrom=? LIMIT 1");
	}

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, iter_sql, -1, &self->iter_stmt, NULL);
		sqlite3_prepare_v2(self->index_db, item_sql, -1, &self->item_stmt, NULL);
		sqlite3_prepare_v2(self->index_db, in_sql, -1, &self->in_stmt, NULL);
	);
	sqlite3_free(iter_sql);
	sqlite3_free(item_sql);
	sqlite3_free(in_sql);
}

void pyfastx_fasta_keys_count(pyfastx_FastaKeys *self) {
	int ret;
	sqlite3_stmt *stmt;

	char *sql = sqlite3_mprintf("SELECT COUNT(1) FROM seq %s %s LIMIT 1",
		self->filter ? "WHERE": "",
		self->filter ? self->filter : ""
	);

	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));
	sqlite3_free(sql);

	PYFASTX_SQLITE_CALL(ret=sqlite3_step(stmt));

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(self->seq_counts = sqlite3_column_int64(stmt, 0));
	} else {
		self->seq_counts = 0;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
}

PyObject *pyfastx_fasta_keys_create(sqlite3 *index_db, uint64_t seq_counts) {
	pyfastx_FastaKeys *keys = PyObject_New(pyfastx_FastaKeys, &pyfastx_FastaKeysType);
	keys->index_db = index_db;
	keys->iter_stmt = NULL;
	keys->item_stmt = NULL;
	keys->in_stmt = NULL;
	keys->seq_counts = seq_counts;
	keys->order = NULL;
	keys->filter = NULL;
	keys->temp_filter = NULL;

	//prepare sql
	pyfastx_fasta_keys_prepare(keys);

	//Py_INCREF(keys);
	return (PyObject *)keys;
}

void pyfastx_fasta_keys_dealloc(pyfastx_FastaKeys *self) {
	if (self->iter_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->iter_stmt));
		self->iter_stmt = NULL;
	}

	if (self->item_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->item_stmt));
		self->item_stmt = NULL;
	}

	if (self->in_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->in_stmt));
		self->in_stmt = NULL;
	}

	if (self->temp_filter) {
		free(self->temp_filter);
		self->temp_filter = NULL;
	}

	if (self->filter) {
		free(self->filter);
		self->filter = NULL;
	}

	if (self->order) {
		sqlite3_free(self->order);
		self->order = NULL;
	}

	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *pyfastx_fasta_keys_iter(pyfastx_FastaKeys *self) {
	pyfastx_fasta_keys_prepare(self);

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fasta_keys_next(pyfastx_FastaKeys *self) {
	int ret;

	PYFASTX_SQLITE_CALL(ret=sqlite3_step(self->iter_stmt));

	if (ret == SQLITE_ROW) {
		char *p;
		int nbytes;
		PYFASTX_SQLITE_CALL(nbytes = sqlite3_column_bytes(self->iter_stmt, 0));
		PyObject *name = PyUnicode_New(nbytes, 255);
		Py_UCS1 *s = PyUnicode_1BYTE_DATA(name);
		PYFASTX_SQLITE_CALL(p = (char *)sqlite3_column_text(self->iter_stmt, 0));
		memcpy(s, p, nbytes);
		return name;
	}

	return NULL;
}

uint64_t pyfastx_fasta_keys_length(pyfastx_FastaKeys *self) {
	return self->seq_counts;
}

PyObject *pyfastx_fasta_keys_repr(pyfastx_FastaKeys *self) {
	return PyUnicode_FromFormat("<FastaKeys> contains %ld keys", self->seq_counts);
}

PyObject *pyfastx_fasta_keys_subscript(pyfastx_FastaKeys *self, PyObject *item) {
	int ret;

	if (PyIndex_Check(item)) {
		Py_ssize_t i;

		i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		if (i == -1 && PyErr_Occurred()) {
			return NULL;
		}

		if (i < 0){
			i = i + self->seq_counts;
		}

		++i;

		if (i > self->seq_counts){
			PyErr_SetString(PyExc_IndexError, "index out of range");
			return NULL;
		}

		if (self->filter || self->order) {
			--i;
		}

		PYFASTX_SQLITE_CALL(
			sqlite3_reset(self->item_stmt);
			sqlite3_bind_int(self->item_stmt, 1, i);
			ret = sqlite3_step(self->item_stmt);
		);

		if (ret == SQLITE_ROW) {
			char *p;
			int nbytes;
			PYFASTX_SQLITE_CALL(nbytes = sqlite3_column_bytes(self->item_stmt, 0));
			PyObject *name = PyUnicode_New(nbytes, 255);
			Py_UCS1 *s = PyUnicode_1BYTE_DATA(name);
			PYFASTX_SQLITE_CALL(p = (char *)sqlite3_column_text(self->item_stmt, 0));
			memcpy(s, p, nbytes);
			return name;
		} else {
			PyErr_Format(PyExc_ValueError, "get item error, code: %d", ret);
			return NULL;
		}
	} else if (PySlice_Check(item)) {
		Py_ssize_t slice_start, slice_stop, slice_step, slice_len;

		if (PySlice_Unpack(item, &slice_start, &slice_stop, &slice_step) < 0) {
			return NULL;
		}
		slice_len = PySlice_AdjustIndices(self->seq_counts, &slice_start, &slice_stop, slice_step);

		if (slice_len <= 0) {
			return PyList_New(0);
		}

		char *sql = sqlite3_mprintf("SELECT chrom FROM seq %s %s %s LIMIT %d OFFSET %d",
			self->filter ? "WHERE" : "",
			self->filter ? self->filter : "",
			self->order ? self->order : "ORDER BY ID",
			slice_len, slice_start
		);

		sqlite3_stmt *stmt;
		char *p;
		PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));
		sqlite3_free(sql);

		PyObject *chroms = PyList_New(0);
		PYFASTX_SQLITE_CALL(ret = sqlite3_step(stmt));
		while (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(p = (char *)sqlite3_column_text(stmt, 0));
			PyObject *chrom = Py_BuildValue("s", p);
			PyList_Append(chroms, chrom);
			Py_DECREF(chrom);
			PYFASTX_SQLITE_CALL(ret = sqlite3_step(stmt));
		}

		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		return chroms;
	} else {
		PyErr_Format(PyExc_TypeError, "fakeys indices must be integers or slices");
		return NULL;
	}
}

int pyfastx_fasta_keys_contains(pyfastx_FastaKeys *self, PyObject *key) {
	char *name;
	int ret;

	if (!PyUnicode_CheckExact(key)) {
		return 0;
	}

	name = (char *)PyUnicode_AsUTF8(key);

	PYFASTX_SQLITE_CALL(
		sqlite3_bind_text(self->in_stmt, 1, name, -1, NULL);
		ret = sqlite3_step(self->in_stmt);
	);

	return ret==SQLITE_ROW ? 1 : 0;
}

PyObject *pyfastx_fasta_keys_sort(pyfastx_FastaKeys *self, PyObject *args, PyObject *kwargs) {
	char *by = "id";
	int reverse = 0;
	int sort = 0;
	
	static char* kwlist[] = {"by", "reverse", NULL};

	// cannot use uint16_t to parse Python bool, should use int declare
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|si", kwlist, &by, &reverse)) {
		return NULL;
	}

	//set sort column
	if (strcmp(by, "id") == 0) {
		sort = 0;
	} else if (strcmp(by, "name") == 0) {
		sort = 1;
	} else if (strcmp(by, "length") == 0) {
		sort = 2;
	} else {
		PyErr_SetString(PyExc_ValueError, "key only can be id, name or length");
		return NULL;
	}

	if (sort || reverse) {
		self->order = sqlite3_mprintf("ORDER BY %s %s", SORTS[sort], ORDERS[reverse]);
	}

	pyfastx_fasta_keys_prepare(self);

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fasta_keys_richcompare(pyfastx_FastaKeys *self, PyObject *other, int op) {
	char *when;
	long slen;
	char signs[][3] = {"<", "<=", "=", "<>", ">", ">="};
	int signt = 0;

	if (!PyLong_Check(other)) {
		PyErr_SetString(PyExc_ValueError, "the compared item must be an integer");
		return NULL;
	}

	slen = PyLong_AsLong(other);
	
	switch (op) {
		case Py_LT: signt = 0; break;
		case Py_LE: signt = 1; break;
		case Py_EQ: signt = 2; break;
		case Py_NE: signt = 3; break;
		case Py_GT: signt = 4; break;
		case Py_GE: signt = 5; break;
	}

	if (self->temp_filter) {
		when = sqlite3_mprintf(" AND slen %s %ld", signs[signt], slen);
		self->temp_filter = (char *)realloc(self->temp_filter, strlen(self->temp_filter) + strlen(when) + 1);
		strcat(self->temp_filter, when);
		sqlite3_free(when);
	} else {
		when = sqlite3_mprintf("slen %s %ld", signs[signt], slen);
		self->temp_filter = (char *)malloc(strlen(when) + 1);
		strcpy(self->temp_filter, when);
		sqlite3_free(when);
	}

	return Py_BuildValue("s", self->temp_filter);
}

PyObject *pyfastx_fasta_keys_like(pyfastx_FastaKeys *self, PyObject *tag) {
	if (!PyUnicode_CheckExact(tag)) {
		PyErr_SetString(PyExc_ValueError, "the tag after % must be a string");
		return NULL;
	}

	return PyUnicode_FromFormat("chrom LIKE '%%%s%%'", (char *)PyUnicode_AsUTF8(tag));
}

PyObject *pyfastx_fasta_keys_filter(pyfastx_FastaKeys *self, PyObject *args) {
	Py_ssize_t l;
	Py_ssize_t c = PyTuple_Size(args);

	if (c <= 0) {
		PyErr_SetString(PyExc_ValueError, "no comparison condition provided");
		return NULL;
	}

	PyObject *sep = Py_BuildValue("s", " AND ");
	PyObject *cat = PyUnicode_Join(sep, args);
	char *tmp = (char *)PyUnicode_AsUTF8AndSize(cat, &l);

	if (self->filter) {
		self->filter = (char *)realloc(self->filter, l+1);
	} else {
		self->filter = (char *)malloc(l+1);
	}
	
	strcpy(self->filter, tmp);
	Py_DECREF(sep);
	Py_DECREF(cat);
	
	if (self->temp_filter) {
		free(self->temp_filter);
		self->temp_filter = NULL;
	}

	pyfastx_fasta_keys_count(self);
	pyfastx_fasta_keys_prepare(self);

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fasta_keys_reset(pyfastx_FastaKeys *self) {
	sqlite3_stmt *stmt;
	int ret;

	if (self->filter) {
		free(self->filter);
		self->filter = NULL;
	}

	if (self->temp_filter) {
		free(self->temp_filter);
		self->temp_filter = NULL;
	}

	if (self->order) {
		sqlite3_free(self->order);
		self->order = NULL;
	}

	pyfastx_fasta_keys_prepare(self);

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, "SELECT seqnum FROM stat", -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			self->seq_counts = sqlite3_column_int64(stmt, 0);
			sqlite3_finalize(stmt);
		);
	} else {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_SetString(PyExc_RuntimeError, "get sequence counts error");
		return NULL;
	}

	Py_INCREF(self);
	return (PyObject *)self;
}

static PyMethodDef pyfastx_fasta_keys_methods[] = {
	{"sort", (PyCFunction)pyfastx_fasta_keys_sort, METH_VARARGS|METH_KEYWORDS, NULL},
	{"filter", (PyCFunction)pyfastx_fasta_keys_filter, METH_VARARGS, NULL},
	{"reset", (PyCFunction)pyfastx_fasta_keys_reset, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

//as a number
static PyNumberMethods fasta_keys_as_number = {
	0,
	0,
	0,
	(binaryfunc)pyfastx_fasta_keys_like,
	0,
};

//as a mapping
static PyMappingMethods fasta_keys_as_mapping = {
	0, /* mp_length */
	(binaryfunc)pyfastx_fasta_keys_subscript, /* mp_subscript */
	0, /* mp_ass_subscript */
};

//as a list
static PySequenceMethods fasta_keys_as_sequence = {
	(lenfunc)pyfastx_fasta_keys_length, /*sq_length*/
	0, /*sq_concat*/
	0, /*sq_repeat*/
	0, /*sq_item*/
	0, /* sq_slice */
	0, /* sq_ass_item*/
	0, /* sq_ass_slice */
	(objobjproc)pyfastx_fasta_keys_contains, /*sq_contains*/
	0, /*sq_inplace_concat*/
	0, /*sq_inplace_repeat*/
};

PyTypeObject pyfastx_FastaKeysType = {
    //PyVarObject_HEAD_INIT(&PyType_Type, 0)
    PyVarObject_HEAD_INIT(NULL, 0)
    "Identifier",                     /* tp_name */
    sizeof(pyfastx_FastaKeys),       /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_fasta_keys_dealloc,                              /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)pyfastx_fasta_keys_repr,                              /* tp_repr */
    &fasta_keys_as_number,                              /* tp_as_number */
    &fasta_keys_as_sequence,                              /* tp_as_sequence */
    &fasta_keys_as_mapping,   /* tp_as_mapping */
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
    (richcmpfunc)pyfastx_fasta_keys_richcompare,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    (getiterfunc)pyfastx_fasta_keys_iter,                              /* tp_iter */
    (iternextfunc)pyfastx_fasta_keys_next,                              /* tp_iternext */
    pyfastx_fasta_keys_methods,       /* tp_methods */
    0,       /* tp_members */
    0,       /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    PyType_GenericNew,           /* tp_new */
};