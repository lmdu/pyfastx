#ifndef PYFASTX_FASTQ_KEYS_H
#define PYFASTX_FASTQ_KEYS_H
#include "Python.h"
#include "stdint.h"
#include "sqlite3.h"
#include "util.h"

typedef struct {
	PyObject_HEAD

	sqlite3* index_db;
	sqlite3_stmt *stmt;

	uint64_t read_counts;
} pyfastx_FastqKeys;

extern PyTypeObject pyfastx_FastqKeysType;

#endif