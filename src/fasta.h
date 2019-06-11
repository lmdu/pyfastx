#ifndef PYFASTX_FASTA_H
#define PYFASTX_FASTA_H
#include <zlib.h>
#include <Python.h>
#include <sqlite3.h>
#include "kseq.h"
#include "zran.h"

KSEQ_DECLARE(gzFile)

//make sequence iterator
typedef struct {
	PyObject_HEAD

	//fasta or fastq file path and name
	char* file_name;

	//sqlite3 index file path and name
	char* index_file;
	
	//gzip open file handle
	gzFile gzfd;

	//file open handle
	FILE* fd;
	
	//kseqs for reading from fasta
	kseq_t* kseqs;

	//always output uppercase sequence
	int uppercase;

	//is gzip compressed file
	//0 not gzip file
	//1 is gzip file
	int gzip_format;

	//total sequence counts
	int seq_counts;

	//total sequence length (bp)
	long seq_length;

	//a float for GC content (%)
	float gc_content;

	//a dict for storing A T G C N (unknown base) counts in fasta
	PyObject *composition;

	//sqlite3 handle for index
	sqlite3* index_db;

	//gzip random access index
	zran_index_t *gzip_index;

} pyfastx_Fasta;

extern PyTypeObject pyfastx_FastaType;

PyObject* test(pyfastx_Fasta *self, PyObject *args);
PyObject* fastx_tp_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);
PyObject* fastx_tp_iter(pyfastx_Fasta *self);
PyObject* fastx_tp_next(pyfastx_Fasta *self);
int fastx_get_item(pyfastx_Fasta *self, PyObject *key);
PyObject *fastx_get_key(pyfastx_Fasta *self, PyObject *key);
void fastx_tp_dealloc(pyfastx_Fasta *self);
int fastx_get_len(pyfastx_Fasta *self);
int fastx_get_val(pyfastx_Fasta *self, PyObject *key, PyObject *val);

void _pyfastx_build_gzip_index(pyfastx_Fasta *self);
void _pyfastx_load_gzip_index(pyfastx_Fasta *self);
PyObject *pyfastx_build_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs);
PyObject *get_sub_seq(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs);

#endif