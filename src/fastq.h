#ifndef PYFASTX_FASTQ_H
#define PYFASTX_FASTQ_H

#include "Python.h"
#include "stdint.h"
#include "zlib.h"
#include "kseq.h"
#include "zran.h"
#include "util.h"
#include "sqlite3.h"

typedef struct {
	PyObject_HEAD

	//fastq file name
	char* file_name;

	//index file
	char* index_file;

	//phred
	uint16_t phred;

	//total read counts
	uint64_t read_counts;

	//total sequence length;
	uint64_t seq_length;

	//GC content
	float gc_content;

	//base composition
	//PyObject* composition;

	//is gzip file
	uint16_t gzip_format;

	//sqlite3 connection to index file
	sqlite3* index_db;

	//file handle for zran index
	FILE* fd;

	//gzip file handle
	gzFile gzfd;

	//kstream for reading line from fastq
	kstream_t *ks;

	//kseq for iteration
	kseq_t *kseq;

	//the ID of the current iter read
	//uint64_t iter_id;
	sqlite3_stmt *iter_stmt;

	//if build_index is True means has index
	uint16_t has_index;

	//gzip index
	zran_index_t* gzip_index;


} pyfastx_Fastq;

extern PyTypeObject pyfastx_FastqType;

void pyfastx_fastq_dealloc(pyfastx_Fastq *self);
void pyfastx_fastq_create_index(pyfastx_Fastq *self);
void pyfastx_fastq_calc_composition(pyfastx_Fastq *self);
PyObject *pyfastx_fastq_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);

#endif