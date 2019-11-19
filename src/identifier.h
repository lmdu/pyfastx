#ifndef PYFASTX_IDENTIFIER_H
#define PYFASTX_IDENTIFIER_H
#include "Python.h"
#include <stdint.h>
#include "sqlite3.h"
#include "pycompat.h"

typedef struct {
	PyObject_HEAD

	//index for fast random access to sequence
	sqlite3* index_db;

	//sqlite3 handle
	sqlite3_stmt *stmt;

	//sequence counts
	uint32_t seq_counts;

	//file format 1: fasta, 2: fastq
	//uint16_t format;

	//sort by 1: id, 2: name, 3: length
	char *sort;

	//oder 0: asc, 1: desc
	char *order;

	//filter string
	char *temp_filter;
	char *filter;

} pyfastx_Identifier;

extern PyTypeObject pyfastx_IdentifierType;

PyObject *pyfastx_identifier_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);
PyObject *pyfastx_identifier_sort(pyfastx_Identifier *self, PyObject *args, PyObject *kwargs);
PyObject *pyfastx_identifier_iter(pyfastx_Identifier *self);
PyObject *pyfastx_identifier_next(pyfastx_Identifier *self);
int pyfastx_identifier_length(pyfastx_Identifier *self);
PyObject *pyfastx_identifier_repr(pyfastx_Identifier *self);
PyObject *pyfastx_identifier_item(pyfastx_Identifier *self, Py_ssize_t i);
int pyfastx_identifier_contains(pyfastx_Identifier *self, PyObject *key);

#endif
