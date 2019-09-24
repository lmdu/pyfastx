#ifndef PYFASTX_IDENTIFIER_H
#define PYFASTX_IDENTIFIER_H
#include "Python.h"
#include "sqlite3.h"

typedef struct {
	PyObject_HEAD

	//index for fast random access to sequence
	sqlite3* index_db;

	//sqlite3 handle
	sqlite3_stmt *stmt;

	//sequence counts
	uint32_t seq_counts;

} pyfastx_Identifier;

extern PyTypeObject pyfastx_IdentifierType;

PyObject *pyfastx_identifier_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);
PyObject *pyfastx_identifier_iter(pyfastx_Identifier *self);
PyObject *pyfastx_identifier_next(pyfastx_Identifier *self);
int pyfastx_identifier_length(pyfastx_Identifier *self);
PyObject *pyfastx_identifier_repr(pyfastx_Identifier *self);
PyObject *pyfastx_identifier_item(pyfastx_Identifier *self, Py_ssize_t i);
int pyfastx_identifier_contains(pyfastx_Identifier *self, PyObject *key);

#endif
