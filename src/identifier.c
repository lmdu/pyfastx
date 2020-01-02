#include "identifier.h"
#include "util.h"

/*PyObject *pyfastx_identifier_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	pyfastx_Identifier *obj = (pyfastx_Identifier *)type->tp_alloc(type, 0);
	return (PyObject *)obj;
}*/

/*
void pyfastx_identifier_dealloc(pyfastx_Identifier *self){
	if (self->stmt != NULL) {
		sqlite3_finalize(self->stmt);
	}

	Py_TYPE(self)->tp_free(self);
}*/

PyObject *pyfastx_identifier_iter(pyfastx_Identifier *self) {
	PyObject *tmp;
	char *sql;

	if (self->filter) {
		tmp = PyUnicode_FromFormat("SELECT chrom FROM seq WHERE %s ORDER BY %s %s",
			self->filter, self->sort, self->order);
	} else {
		tmp = PyUnicode_FromFormat("SELECT chrom FROM seq ORDER BY %s %s",
			self->sort, self->order);
	}

	sql = PyUnicode_AsUTF8(tmp);

	sqlite3_prepare_v2(self->index_db, sql, -1, &self->stmt, NULL);

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_identifier_next(pyfastx_Identifier *self){
	int nbytes;
	char *name;

	if (sqlite3_step(self->stmt) != SQLITE_ROW){
		sqlite3_reset(self->stmt);
		return NULL;
	}

	nbytes = sqlite3_column_bytes(self->stmt, 0);
	name = (char *)malloc(nbytes + 1);
	memcpy(name, (char *)sqlite3_column_text(self->stmt, 0), nbytes);
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
	PyObject *tmp;
	char *sql;
	
	if(i < 0){
		i = i + self->seq_counts;
	}

	if(i >= self->seq_counts){
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return NULL;
	}

	if (self->filter) {
		tmp = PyUnicode_FromFormat("%s %s ORDER BY %s %s LIMIT 1 OFFSET %ld",
			"SELECT chrom FROM seq WHERE", self->filter, self->sort, self->order, i);
	} else {
		tmp = PyUnicode_FromFormat("%s ORDER BY %s %s LIMIT 1 OFFSET %ld",
			"SELECT chrom FROM seq", self->sort, self->order, i);
	}

	sql = PyUnicode_AsUTF8(tmp);
	
	sqlite3_prepare_v2(self->index_db, sql, -1, &self->stmt, NULL);

	if (sqlite3_step(self->stmt) == SQLITE_ROW) {
		int nbytes = sqlite3_column_bytes(self->stmt, 0);
		char *name = (char *)malloc(nbytes + 1);
		memcpy(name, (char *)sqlite3_column_text(self->stmt, 0), nbytes);
		name[nbytes] = '\0';
		sqlite3_finalize(self->stmt);
		return Py_BuildValue("s", name);
	} else {
		PyErr_SetString(PyExc_ValueError, "get item error");
		return NULL;
	}
}

int pyfastx_identifier_contains(pyfastx_Identifier *self, PyObject *key){
	char *name;
	const char *sql;

	if(!PyUnicode_CheckExact(key)){
		return 0;
	}

	name = PyUnicode_AsUTF8(key);

	sql = "SELECT * FROM seq WHERE chrom=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &self->stmt, NULL);
	sqlite3_bind_text(self->stmt, 1, name, -1, NULL);
	if(sqlite3_step(self->stmt) != SQLITE_ROW){
		sqlite3_reset(self->stmt);
		return 0;
	} else {
		sqlite3_reset(self->stmt);
		return 1;
	}
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
		self->sort = "ID";
	} else if (strcmp(key, "name") == 0) {
		self->sort = "chrom";
	} else if (strcmp(key, "length") == 0) {
		self->sort = "slen";
	} else {
		PyErr_SetString(PyExc_ValueError, "key only can be id, name or length");
		return NULL;
	}

	//set sort order
	if (reverse) {
		self->order = "DESC";
	} else {
		self->order = "ASC";
	}

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_idnetifier_richcompare(pyfastx_Identifier *self, PyObject *other, int op) {
	char when[100];
	char *sign;
	uint32_t slen;
	int flen;

	if (!PyLong_Check(other)) {
		PyErr_SetString(PyExc_ValueError, "the compared item must be an integer");
		return NULL;
	}

	memset(when, '\0', sizeof(when));

	slen = PyLong_AsLong(other);
	
	switch (op) {
		case Py_LT: sign = "<"; break;
		case Py_LE: sign = "<="; break;
		case Py_EQ: sign = "="; break;
		case Py_NE: sign = "<>"; break;
		case Py_GT: sign = ">"; break;
		case Py_GE: sign = ">="; break;
		default:
			PyErr_SetString(PyExc_ValueError, "the comparison symbol not supported");
			return Py_NotImplemented;
	}

	if (self->temp_filter == NULL) {
		flen = sprintf(when, "slen %s %d", sign, slen);
		self->temp_filter = (char *)malloc(flen + 1);
		strcpy(self->temp_filter, when);
	} else {
		flen = sprintf(when, " AND slen %s %d", sign, slen);
		self->temp_filter = (char *)realloc(self->temp_filter, strlen(self->temp_filter) + flen + 1);
		strcat(self->temp_filter, when);
	}

	return Py_BuildValue("s", self->temp_filter);
}

PyObject *pyfastx_identifier_like(pyfastx_Identifier *self, PyObject *tag) {
	char *name;

	if (!PyUnicode_CheckExact(tag)) {
		PyErr_SetString(PyExc_ValueError, "the tag after % must be a string");
		return NULL;
	}
	
	name = PyUnicode_AsUTF8(tag);
	
	return PyUnicode_FromFormat("chrom LIKE '%%%s%%'", name);
}

PyObject *pyfastx_identifier_filter(pyfastx_Identifier *self, PyObject *args) {
	PyObject *res;
	PyObject *tmp;
	char *sql;

	Py_ssize_t c = PyTuple_Size(args);

	if (c == 0) {
		PyErr_SetString(PyExc_ValueError, "no comparison condition provided");
		return NULL;
	}

	res = PyUnicode_Join(Py_BuildValue("s", " AND "), args);
	self->filter = PyUnicode_AsUTF8(res);
	
	if (self->temp_filter) {
		free(self->temp_filter);
		self->temp_filter = NULL;
	}

	tmp = PyUnicode_FromFormat("SELECT COUNT(*) FROM seq WHERE %s", self->filter);
	sql = PyUnicode_AsUTF8(tmp);

	sqlite3_prepare_v2(self->index_db, sql, -1, &self->stmt, NULL);

	if (sqlite3_step(self->stmt) == SQLITE_ROW) {
		self->seq_counts = sqlite3_column_int(self->stmt, 0);
		sqlite3_finalize(self->stmt);
	} else {
		self->seq_counts = 0;
	}

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_identifier_reset(pyfastx_Identifier *self) {
	self->sort = NULL;
	self->order = NULL;
	self->filter = NULL;

	if (self->temp_filter) {
		free(self->temp_filter);
		self->temp_filter = NULL;
	}
	
	sqlite3_prepare_v2(self->index_db, "SELECT COUNT(*) FROM seq", -1, &self->stmt, NULL);

	if (sqlite3_step(self->stmt) == SQLITE_ROW) {
		self->seq_counts = sqlite3_column_int(self->stmt, 0);
		sqlite3_finalize(self->stmt);
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
    0,                              /* tp_dealloc */
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
    0,           /* tp_new */
};