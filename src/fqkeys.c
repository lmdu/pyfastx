#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "fqkeys.h"

PyObject *pyfastx_fastq_keys_create(sqlite3 *index_db, Py_ssize_t read_counts) {
	pyfastx_FastqKeys *keys = PyObject_New(pyfastx_FastqKeys, &pyfastx_FastqKeysType);
	keys->index_db = index_db;
	keys->read_counts = read_counts;
	keys->iter_stmt = NULL;
	keys->item_stmt = NULL;
	keys->in_stmt = NULL;

	//prepare sql
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(keys->index_db, "SELECT name FROM read ORDER BY ID", -1, &keys->iter_stmt, NULL);
		sqlite3_prepare_v2(keys->index_db, "SELECT name FROM read WHERE ID=? LIMIT 1", -1, &keys->item_stmt, NULL);
		sqlite3_prepare_v2(keys->index_db, "SELECT 1 FROM read WHERE name=? LIMIT 1", -1, &keys->in_stmt, NULL);
	);

	return (PyObject *)keys;
}

void pyfastx_fastq_keys_dealloc(pyfastx_FastqKeys *self) {
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

	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *pyfastx_fastq_keys_repr(pyfastx_FastqKeys *self) {
	return PyUnicode_FromFormat("<FastqKeys> contains %ld keys", self->read_counts);
}

PyObject *pyfastx_fastq_keys_iter(pyfastx_FastqKeys *self) {
	if (self->iter_stmt) {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(self->iter_stmt));
		self->iter_stmt = NULL;
	}
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index_db, "SELECT name FROM read ORDER BY ID", -1, &self->iter_stmt, NULL);
	);
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fastq_keys_next(pyfastx_FastqKeys *self) {
	int ret;
	char *name;

	PYFASTX_SQLITE_CALL(ret=sqlite3_step(self->iter_stmt));

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(name = (char *)sqlite3_column_text(self->iter_stmt, 0));
		return Py_BuildValue("s", name);
	}

	return NULL;
}

Py_ssize_t pyfastx_fastq_keys_length(pyfastx_FastqKeys *self) {
	return self->read_counts;
}

PyObject *pyfastx_fastq_keys_item(pyfastx_FastqKeys *self, Py_ssize_t i) {
	int ret;
	char *name;

	if (i < 0){
		i = i + self->read_counts;
	}

	++i;

	if (i > self->read_counts){
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return NULL;
	}

	PYFASTX_SQLITE_CALL(
		sqlite3_reset(self->item_stmt);
		sqlite3_bind_int64(self->item_stmt, 1, i);
		ret = sqlite3_step(self->item_stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(name = (char *)sqlite3_column_text(self->item_stmt, 0));
		return Py_BuildValue("s", name);
	} else {
		PyErr_Format(PyExc_ValueError, "get item error, code: %d", ret);
		return NULL;
	}
}

int pyfastx_fastq_keys_contains(pyfastx_FastqKeys *self, PyObject *key) {
	int ret;
	char *name;

	if (!PyUnicode_CheckExact(key)) {
		return 0;
	}

	name = (char *)PyUnicode_AsUTF8(key);

	PYFASTX_SQLITE_CALL(
		sqlite3_reset(self->in_stmt);
		sqlite3_bind_text(self->in_stmt, 1, name, -1, NULL);
		ret = sqlite3_step(self->in_stmt);
	);

	return ret==SQLITE_ROW ? 1 : 0;
}

//as a list
static PySequenceMethods fastq_keys_as_sequence = {
	.sq_length = (lenfunc)pyfastx_fastq_keys_length,
	.sq_item = (ssizeargfunc)pyfastx_fastq_keys_item,
	.sq_contains = (objobjproc)pyfastx_fastq_keys_contains,
};

PyTypeObject pyfastx_FastqKeysType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "FastqKeys",
    .tp_basicsize = sizeof(pyfastx_FastqKeys),
    .tp_dealloc = (destructor)pyfastx_fastq_keys_dealloc,
    .tp_repr = (reprfunc)pyfastx_fastq_keys_repr,
    .tp_as_sequence = &fastq_keys_as_sequence,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_iter = (getiterfunc)pyfastx_fastq_keys_iter,
    .tp_iternext = (iternextfunc)pyfastx_fastq_keys_next,
    .tp_alloc = PyType_GenericAlloc,
    .tp_new = PyType_GenericNew,
};