#ifndef PYFASTX_FASTA_KEYS_H
#define PYFASTX_FASTA_KEYS_H
#include "Python.h"
#include <stdint.h>
#include "sqlite3.h"

typedef struct {
	PyObject_HEAD

	//index for fast random access to sequence
	sqlite3* index_db;

	//sqlite3 handle
	sqlite3_stmt *stmt;

	//sequence counts
	uint64_t seq_counts;

	//file format 1: fasta, 2: fastq
	//uint16_t format;

	//sort by 0: id, 1: name, 2: length
	uint16_t sort;

	//oder 0: asc, 1: desc
	uint16_t order;

	//need to update results
	uint8_t update;

	//filter string
	char *temp_filter;
	char *filter;

} pyfastx_FastaKeys;

extern PyTypeObject pyfastx_FastaKeysType;

#endif
