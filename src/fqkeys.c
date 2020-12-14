#include "fqkeys.h"

void pyfastx_fastq_keys_dealloc(pyfastx_FastqKeys *self) {
	if (self->stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->stmt));
	}

	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *pyfastx_fastq_keys_repr(pyfastx_FastqKeys *self) {
	return PyUnicode_FromFormat("<FastqKeys> contains %ld keys", self->read_counts);
}

PyObject *pyfastx_fastq_keys_iter(pyfastx_FastqKeys *self) {
	if (self->stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->stmt));
		self->stmt = NULL;
	}
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, "SELECT name FROM read ORDER BY ID", -1, &self->stmt, NULL));
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fastq_keys_next(pyfastx_FastqKeys *self){
	int nbytes;
	char *name;
	int ret;

	PYFASTX_SQLITE_CALL(ret=sqlite3_step(self->stmt));

	if (ret == SQLITE_ROW){
		PYFASTX_SQLITE_CALL(
			nbytes = sqlite3_column_bytes(self->stmt, 0);
			name = (char *)malloc(nbytes + 1);
			memcpy(name, (char *)sqlite3_column_text(self->stmt, 0), nbytes);
		);
		name[nbytes] = '\0';
		return Py_BuildValue("s", name);
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(self->stmt));
	self->stmt = NULL;
	return NULL;
}

uint64_t pyfastx_fastq_keys_length(pyfastx_FastqKeys *self) {
	return self->read_counts;
}

PyObject *pyfastx_fastq_keys_item(pyfastx_FastqKeys *self, Py_ssize_t i) {
	sqlite3_stmt *stmt;
	char *sql;
	int ret;
	int nbytes;
	char *name;

	if (i < 0){
		i = i + self->read_counts;
	}

	++i;

	if (i > self->read_counts){
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return NULL;
	}

	sql = "SELECT name FROM read WHERE ID=?";
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL));
	PYFASTX_SQLITE_CALL(
		sqlite3_bind_int64(stmt, 1, i);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(nbytes = sqlite3_column_bytes(stmt, 0));
		name = (char *)malloc(nbytes + 1);
		PYFASTX_SQLITE_CALL(memcpy(name, (char *)sqlite3_column_text(stmt, 0), nbytes));
		name[nbytes] = '\0';
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		return Py_BuildValue("s", name);
	} else {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_Format(PyExc_ValueError, "get item error, code: %d", ret);
		return NULL;
	}
}

int pyfastx_fastq_keys_contains(pyfastx_FastqKeys *self, PyObject *key) {
	char *name;
	sqlite3_stmt *stmt;
	char *sql;
	int ret;

	if (!PyUnicode_CheckExact(key)) {
		return 0;
	}

	name = (char *)PyUnicode_AsUTF8(key);
	sql = sqlite3_mprintf("SELECT 1 FROM read WHERE name=? LIMIT 1");

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	return ret==SQLITE_ROW ? 1 : 0;
}

//as a list
static PySequenceMethods fastq_keys_as_sequence = {
	(lenfunc)pyfastx_fastq_keys_length, /*sq_length*/
	0, /*sq_concat*/
	0, /*sq_repeat*/
	(ssizeargfunc)pyfastx_fastq_keys_item, /*sq_item*/
	0, /*sq_slice */
	0, /*sq_ass_item*/
	0, /*sq_ass_splice*/
	(objobjproc)pyfastx_fastq_keys_contains, /*sq_contains*/
	0, /*sq_inplace_concat*/
	0, /*sq_inplace_repeat*/
};

PyTypeObject pyfastx_FastqKeysType = {
    //PyVarObject_HEAD_INIT(&PyType_Type, 0)
    PyVarObject_HEAD_INIT(NULL, 0)
    "FastqKeys",                     /* tp_name */
    sizeof(pyfastx_FastqKeys),       /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_fastq_keys_dealloc,                              /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)pyfastx_fastq_keys_repr,                              /* tp_repr */
    0,                              /* tp_as_number */
    &fastq_keys_as_sequence,                              /* tp_as_sequence */
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
    (getiterfunc)pyfastx_fastq_keys_iter,                              /* tp_iter */
    (iternextfunc)pyfastx_fastq_keys_next,                              /* tp_iternext */
    0,       /* tp_methods */
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