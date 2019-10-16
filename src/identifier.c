#include "identifier.h"

PyObject *pyfastx_identifier_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	pyfastx_Identifier *obj = (pyfastx_Identifier *)type->tp_alloc(type, 0);
	return (PyObject *)obj;
}

void pyfastx_identifier_dealloc(pyfastx_Identifier *self){
	if(self->stmt != NULL){
		sqlite3_finalize(self->stmt);
	}
}

PyObject *pyfastx_identifier_iter(pyfastx_Identifier *self){
	char *key;
	char *order;
	char *sql;

	if (self->sort == 1) {
		key = "ID";
	} else if (self->sort == 2) {
		key = "chrom";
	} else if (self->sort == 3) {
		key = "slen";
	} else {
		key = NULL;
	}

	if (self->order == 0) {
		order = "ASC";
	} else {
		order = "DESC";
	}

	sql = (char *)malloc(50 * sizeof(char));
	sprintf(sql, "SELECT chrom FROM seq ORDER BY %s %s;", key, order);

	sqlite3_prepare_v2(self->index_db, sql, -1, &self->stmt, NULL);
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_identifier_next(pyfastx_Identifier *self){
	if(sqlite3_step(self->stmt) != SQLITE_ROW){
		sqlite3_reset(self->stmt);
		return NULL;
	}
	char *name;
	name = (char *)sqlite3_column_text(self->stmt, 0);
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

	char *name = (char *)sqlite3_column_text(self->stmt, 0);

	sqlite3_reset(self->stmt);

	return Py_BuildValue("s", name);
}

int pyfastx_identifier_contains(pyfastx_Identifier *self, PyObject *key){
	if(!PyUnicode_Check(key)){
		return 0;
	}

	char *name = PyUnicode_AsUTF8(key);

	sqlite3_prepare_v2(self->index_db, "SELECT * FROM seq WHERE chrom=? LIMIT 1;", -1, &self->stmt, NULL);
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
		self->sort = 1;
	} else if (strcmp(key, "name") == 0) {
		self->sort = 2;
	} else if (strcmp(key, "length") == 0) {
		self->sort = 3;
	} else {
		PyErr_SetString(PyExc_ValueError, "key only can be id, name or length");
		return NULL;
	}

	//set sort order
	self->order = reverse;

	return (PyObject *)self;
}

static PyMethodDef pyfastx_identifier_methods[] = {
	{"sort", (PyCFunction)pyfastx_identifier_sort, METH_VARARGS|METH_KEYWORDS},
	{NULL, NULL, 0, NULL}
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
    0,                              /* tp_as_number */
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
    0,                              /* tp_richcompare */
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