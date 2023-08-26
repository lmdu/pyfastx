#ifndef PYFASTX_FASTQ_H
#define PYFASTX_FASTQ_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "zlib.h"
#include "kseq.h"
#include "zran.h"
#include "util.h"
#include "sqlite3.h"

typedef struct {
	PyObject_HEAD

	//phred
	int phred;

	//is gzip file
	int gzip_format;

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
	Py_ssize_t cache_soff;

	//cache end offset
	Py_ssize_t cache_eoff;

	//iteration mode
	int iterating;

	PyObject *fastq;

} pyfastx_FastqMiddleware;


typedef struct {
	PyObject_HEAD

	//fastq file name
	PyObject* file_obj;

	//index file
	char* index_file;

	//total read counts
	Py_ssize_t read_counts;

	//total sequence length;
	Py_ssize_t seq_length;

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
	int has_index;

	//average length
	double avg_length;

	//min max length
	Py_ssize_t maxlen;
	Py_ssize_t minlen;

	//min and max quality score
	int maxqual;
	int minqual;

	//iterate with full name
	int full_name;

	pyfastx_FastqMiddleware* middle;

	PyObject* (*func) (pyfastx_FastqMiddleware *);

} pyfastx_Fastq;

extern PyTypeObject pyfastx_FastqType;

void pyfastx_fastq_calc_composition(pyfastx_Fastq *self);

#endif