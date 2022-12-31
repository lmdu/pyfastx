#ifndef PYFASTX_FASTA_KEYS_H
#define PYFASTX_FASTA_KEYS_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "sqlite3.h"

typedef struct {
	PyObject_HEAD

	//index for fast random access to sequence
	sqlite3* index_db;

	//sqlite3 handle
	sqlite3_stmt *iter_stmt;
	sqlite3_stmt *item_stmt;
	sqlite3_stmt *in_stmt;

	//sequence counts
	Py_ssize_t seq_counts;

	//filter string
	char *temp_filter;
	char *filter;

	//order by
	char *order;

} pyfastx_FastaKeys;

extern PyTypeObject pyfastx_FastaKeysType;

PyObject *pyfastx_fasta_keys_create(sqlite3 *index_db, Py_ssize_t seq_counts);

#endif
