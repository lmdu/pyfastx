#include "identifier.h"
#include "util.h"

PyObject *pyfastx_identifier_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	pyfastx_Identifier *obj = (pyfastx_Identifier *)type->tp_alloc(type, 0);
	return (PyObject *)obj;
}

PyObject *make_new_identifier(pyfastx_Identifier *self, char *filter) {
	pyfastx_Identifier *obj = PyObject_New(pyfastx_Identifier, &pyfastx_IdentifierType);
	
	if (!obj) {
		return NULL;
	}

	obj->index_db = self->index_db;
	obj->filter = filter;
	obj->stmt = self->stmt;
	obj->sort = self->sort;
	obj->order = self->order;

	if (filter) {
		char sql[200] = "SELECT COUNT(*) FROM seq WHERE ";
		if (self->)

		
	} else {
		obj->seq_counts = seq_counts;
	}


	Py_INCREF(obj);
	return (PyObject *)obj;
}

/*
void pyfastx_identifier_dealloc(pyfastx_Identifier *self){
	if (self->stmt != NULL) {
		sqlite3_finalize(self->stmt);
	}

	Py_TYPE(self)->tp_free(self);
}*/

PyObject *pyfastx_identifier_iter(pyfastx_Identifier *self){
	char *key;
	char *order;

	char sql[50];
	sprintf(sql, "SELECT chrom FROM seq ORDER BY %s %s;", self->key, self->order);
	sqlite3_prepare_v2(self->index_db, sql, -1, &self->stmt, NULL);

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_identifier_next(pyfastx_Identifier *self){
	if (sqlite3_step(self->stmt) != SQLITE_ROW){
		sqlite3_reset(self->stmt);
		return NULL;
	}
	int nbytes = sqlite3_column_bytes(self->stmt, 0);
	char *name = (char *)malloc(nbytes + 1);
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
	if(i < 0){
		i = i + self->seq_counts;
	}

	if(i >= self->seq_counts){
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return NULL;
	}

	sqlite3_prepare_v2(self->index_db, "SELECT chrom FROM seq WHERE ID=? LIMIT 1;", -1, &self->stmt, NULL);
	sqlite3_bind_int(self->stmt, 1, i+1);
	sqlite3_step(self->stmt);
	int nbytes = sqlite3_column_bytes(self->stmt, 0);
	char *name = (char *)malloc(nbytes + 1);
	memcpy(name, (char *)sqlite3_column_text(self->stmt, 0), nbytes);
	name[nbytes] = '\0';
	sqlite3_finalize(self->stmt);

	return Py_BuildValue("s", name);
}

int pyfastx_identifier_contains(pyfastx_Identifier *self, PyObject *key){
	if(!PyString_CheckExact(key)){
		return 0;
	}

	char *name = PyUnicode_AsUTF8(key);

	const char *sql = "SELECT * FROM seq WHERE chrom=? LIMIT 1;";
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
		self->sort = "id"
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

	return (PyObject *)self;
}

PyObject *pyfastx_idnetifier_richcompare(PyObject *self, PyObject *other, int op) {
	int64_t l = integer_to_long(other);
	printf("value: %ld\n", l);
	printf("operator: %d\n", op);
	return other;
}

PyObject *pyfastx_identifier_like(PyObject *self, PyObject *tag) {
	if (!PyString_CheckExact(tag)) {
		PyErr_SetString(PyExc_ValueError, "the tag after % must be a string");
		return NULL;
	}

	char *name = PyUnicode_AsUTF8(tag);

	printf("the name is : %s\n", name);
	
	return tag;
}

static PyMethodDef pyfastx_identifier_methods[] = {
	{"sort", (PyCFunction)pyfastx_identifier_sort, METH_VARARGS|METH_KEYWORDS},
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
    pyfastx_identifier_new,           /* tp_new */
};