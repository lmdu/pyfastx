#ifndef PYFASTX_FASTQ_KEYS_H
#define PYFASTX_FASTQ_KEYS_H
#include "Python.h"
#include "stdint.h"
#include "sqlite3.h"
#include "util.h"

typedef struct {
	PyObject_HEAD

	sqlite3* index_db;
	sqlite3_stmt *iter_stmt;
	sqlite3_stmt *item_stmt;
	sqlite3_stmt *in_stmt;

	uint64_t read_counts;
} pyfastx_FastqKeys;

extern PyTypeObject pyfastx_FastqKeysType;

PyObject *pyfastx_fastq_keys_create(sqlite3 *index_db, uint64_t read_counts);

#endif