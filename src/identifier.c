#include "identifier.h"
#include "util.h"

char SORTS[][6] = {"ID", "chrom", "slen"};
char ORDERS[][5] = {"ASC", "DESC"};

void pyfastx_identifier_dealloc(pyfastx_Identifier *self){
	if (self->stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->stmt));
	}

	if (self->temp_filter) {
		free(self->temp_filter);
	}

	if (self->filter) {
		free(self->filter);
	}

	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *pyfastx_identifier_iter(pyfastx_Identifier *self) {
	char *sql;

	if (self->filter) {
		sql = sqlite3_mprintf("SELECT chrom FROM seq WHERE %s ORDER BY %s %s",
			self->filter, SORTS[self->sort], ORDERS[self->order]);
	} else {
		sql = sqlite3_mprintf("SELECT chrom FROM seq ORDER BY %s %s",
			SORTS[self->sort], ORDERS[self->order]);
	}

	if (self->stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->stmt));
		self->stmt = NULL;
	}

	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &self->stmt, NULL));
	sqlite3_free(sql);

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_identifier_next(pyfastx_Identifier *self){
	int nbytes;
	char *name;
	int ret;

	PYFASTX_SQLITE_CALL(ret=sqlite3_step(self->stmt));

	if (ret != SQLITE_ROW){
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->stmt));
		self->stmt = NULL;
		return NULL;
	}

	PYFASTX_SQLITE_CALL(
		nbytes = sqlite3_column_bytes(self->stmt, 0);
		name = (char *)malloc(nbytes + 1);
		memcpy(name, (char *)sqlite3_column_text(self->stmt, 0), nbytes);
	);
	name[nbytes] = '\0';
	return Py_BuildValue("s", name);
}

int pyfastx_identifier_length(pyfastx_Identifier *self){
	return self->seq_counts;
}

PyObject *pyfastx_identifier_repr(pyfastx_Identifier *self){
	return PyUnicode_FromFormat("<Identifier> contains %d identifiers", self->seq_counts);
}

PyObject *pyfastx_identifier_item(pyfastx_Identifier *self, Py_ssize_t i){
	sqlite3_stmt *stmt;
	char *sql;
	int ret;
	int nbytes;
	char *name;
	
	if (i < 0){
		i = i + self->seq_counts;
	}

	if (i >= self->seq_counts){
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return NULL;
	}

	if (self->filter) {
		sql = sqlite3_mprintf("SELECT chrom FROM seq WHERE %s ORDER BY %s %s LIMIT 1 OFFSET %ld",
			self->filter, SORTS[self->sort], ORDERS[self->order], i);
	} else {
		sql = sqlite3_mprintf("SELECT chrom FROM seq ORDER BY %s %s LIMIT 1 OFFSET %ld",
			SORTS[self->sort], ORDERS[self->order], i);
	}
	
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));
	sqlite3_free(sql);

	PYFASTX_SQLITE_CALL(ret = sqlite3_step(stmt));

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(nbytes = sqlite3_column_bytes(stmt, 0));
		name = (char *)malloc(nbytes + 1);
		PYFASTX_SQLITE_CALL(memcpy(name, (char *)sqlite3_column_text(stmt, 0), nbytes));
		name[nbytes] = '\0';
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		return Py_BuildValue("s", name);
	} else {
		PyErr_SetString(PyExc_ValueError, "get item error");
		return NULL;
	}
}

int pyfastx_identifier_contains(pyfastx_Identifier *self, PyObject *key){
	char *name;
	sqlite3_stmt *stmt;
	const char *sql;
	int ret;

	if(!PyUnicode_CheckExact(key)){
		return 0;
	}

	name = PyUnicode_AsUTF8(key);

	sql = "SELECT * FROM seq WHERE chrom=? LIMIT 1;";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	return ret==SQLITE_ROW ? 1 : 0;
}

PyObject *pyfastx_identifier_sort(pyfastx_Identifier *self, PyObject *args, PyObject *kwargs) {
	char *key = "id";
	int reverse = 0;
	
	static char* kwlist[] = {"key", "reverse", NULL};

	// cannot use uint16_t to parse Python bool, should use int declare
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|si", kwlist, &key, &reverse)) {
		return NULL;
	}

	//set sort column
	if (strcmp(key, "id") == 0) {
		self->sort = 0;
	} else if (strcmp(key, "name") == 0) {
		self->sort = 1;
	} else if (strcmp(key, "length") == 0) {
		self->sort = 2;
	} else {
		PyErr_SetString(PyExc_ValueError, "key only can be id, name or length");
		return NULL;
	}

	//set sort order
	self->order = reverse;

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_idnetifier_richcompare(pyfastx_Identifier *self, PyObject *other, int op) {
	char *when;
	uint32_t slen;
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
		when = sqlite3_mprintf(" AND slen %s %d", signs[signt], slen);
		self->temp_filter = (char *)realloc(self->temp_filter, strlen(self->temp_filter) + strlen(when) + 1);
		strcat(self->temp_filter, when);
		sqlite3_free(when);
	} else {
		when = sqlite3_mprintf("slen %s %d", signs[signt], slen);
		self->temp_filter = (char *)malloc(strlen(when) + 1);
		strcpy(self->temp_filter, when);
		sqlite3_free(when);
	}

	return Py_BuildValue("s", self->temp_filter);
}

PyObject *pyfastx_identifier_like(pyfastx_Identifier *self, PyObject *tag) {
	if (!PyUnicode_CheckExact(tag)) {
		PyErr_SetString(PyExc_ValueError, "the tag after % must be a string");
		return NULL;
	}

	return PyUnicode_FromFormat("chrom LIKE '%%%s%%'", PyUnicode_AsUTF8(tag));
}

PyObject *pyfastx_identifier_filter(pyfastx_Identifier *self, PyObject *args) {
	sqlite3_stmt *stmt;
	char *sql;
	char *tmp;
	int ret;
	Py_ssize_t l;
	Py_ssize_t c = PyTuple_Size(args);

	if (!c) {
		PyErr_SetString(PyExc_ValueError, "no comparison condition provided");
		return NULL;
	}

	PyObject *sep = Py_BuildValue("s", " AND ");
	PyObject *cat = PyUnicode_Join(sep, args);
	tmp = PyUnicode_AsUTF8AndSize(cat, &l);

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

	sql = sqlite3_mprintf("SELECT COUNT(*) FROM seq WHERE %s", self->filter);

	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));
	sqlite3_free(sql);

	PYFASTX_SQLITE_CALL(ret=sqlite3_step(stmt));

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(self->seq_counts = sqlite3_column_int(stmt, 0));
	} else {
		self->seq_counts = 0;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_identifier_reset(pyfastx_Identifier *self) {
	sqlite3_stmt *stmt;
	int ret;
	
	self->sort = 0;
	self->order = 0;

	if (self->filter) {
		free(self->filter);
		self->filter = NULL;
	}

	if (self->temp_filter) {
		free(self->temp_filter);
		self->temp_filter = NULL;
	}
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, "SELECT seqnum FROM stat", -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			self->seq_counts = sqlite3_column_int(stmt, 0);
			sqlite3_finalize(stmt);
		);
	} else {
		PyErr_SetString(PyExc_RuntimeError, "get sequence counts error");
		return NULL;
	}

	Py_INCREF(self);
	return (PyObject *)self;
}

static PyMethodDef pyfastx_identifier_methods[] = {
	{"sort", (PyCFunction)pyfastx_identifier_sort, METH_VARARGS|METH_KEYWORDS, NULL},
	{"filter", (PyCFunction)pyfastx_identifier_filter, METH_VARARGS, NULL},
	{"reset", (PyCFunction)pyfastx_identifier_reset, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

//as a number
static PyNumberMethods identifier_as_number = {
	0,
	0,
	0,
	(binaryfunc)pyfastx_identifier_like,
	0,
};

//as a list
static PySequenceMethods identifier_as_sequence = {
	(lenfunc)pyfastx_identifier_length, /*sq_length*/
	0, /*sq_concat*/
	0, /*sq_repeat*/
	(ssizeargfunc)pyfastx_identifier_item, /*sq_item*/
	0, /*sq_slice */
	0, /*sq_ass_item*/
	0, /*sq_ass_splice*/
	(objobjproc)pyfastx_identifier_contains, /*sq_contains*/
	0, /*sq_inplace_concat*/
	0, /*sq_inplace_repeat*/
};

PyTypeObject pyfastx_IdentifierType = {
    //PyVarObject_HEAD_INIT(&PyType_Type, 0)
    PyVarObject_HEAD_INIT(NULL, 0)
    "Identifier",                     /* tp_name */
    sizeof(pyfastx_Identifier),       /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_identifier_dealloc,                              /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)pyfastx_identifier_repr,                              /* tp_repr */
    &identifier_as_number,                              /* tp_as_number */
    &identifier_as_sequence,                              /* tp_as_sequence */
    0,   /* tp_as_mapping */
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
    (richcmpfunc)pyfastx_idnetifier_richcompare,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    (getiterfunc)pyfastx_identifier_iter,                              /* tp_iter */
    (iternextfunc)pyfastx_identifier_next,                              /* tp_iternext */
    pyfastx_identifier_methods,       /* tp_methods */
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