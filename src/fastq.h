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

	//phred
	uint16_t phred;

	//is gzip file
	uint8_t gzip_format;

	//file handle for zran index
	FILE* fd;

	//gzip file handle
	gzFile gzfd;

	//gzip index
	zran_index_t* gzip_index;

	//iteration stmt
	sqlite3_stmt *iter_stmt;

	//kseq for iteration
	kseq_t *kseq;

	//buffer cache
	char* cache_buff;

	//cache start offset
	int64_t cache_soff;

	//cache end offset
	int64_t cache_eoff;

	//iteration mode
	uint8_t iterating;

} pyfastx_FastqMiddleware;

typedef struct {
	PyObject_HEAD

	//fastq file name
	char* file_name;

	//index file
	char* index_file;

	//total read counts
	uint64_t read_counts;

	//total sequence length;
	uint64_t seq_length;

	//GC content
	float gc_content;

	//sqlite3 connection to index file
	sqlite3* index_db;

	//kstream for reading line from fastq
	kstream_t *ks;

	//the ID of the current iter read
	//uint64_t iter_id;
	sqlite3_stmt *id_stmt;
	sqlite3_stmt *name_stmt;

	//if build_index is True means has index
	uint8_t has_index;

	//average length
	float avg_length;

	//min max length
	uint32_t maxlen;
	uint32_t minlen;

	//min and max quality score
	int maxqual;
	int minqual;

	//iterate with full name
	uint8_t full_name;

	pyfastx_FastqMiddleware* middle;

	PyObject* (*func) (pyfastx_FastqMiddleware *);

} pyfastx_Fastq;

extern PyTypeObject pyfastx_FastqType;

void pyfastx_fastq_dealloc(pyfastx_Fastq *self);
void pyfastx_fastq_create_index(pyfastx_Fastq *self);
void pyfastx_fastq_calc_composition(pyfastx_Fastq *self);
PyObject *pyfastx_fastq_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);

#endif