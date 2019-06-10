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
	
	//gzip open file handle
	gzFile gzfp;

	//file open handle
	FILE* fd;
	
	//kseqs for reading from fasta
	kseq_t* kseqs;
	
	//sqlite3 handle for index
	sqlite3* db;
	
	//fasta or fastq file path and name
	char* file_name;

	//sqlite3 index file path and name
	char* index_file;

	//always output uppercase sequence
	int uppercase;

	//is gzip compressed file
	//0 not gzip file
	//1 is gzip file
	int is_gzip;

	//total sequence counts
	int seq_counts;

	//total sequence length (bp)
	int seq_length;

	//GC content (%)
	float gc_content;

	//A T G C N (unknown base) counts in fasta
	int a_counts;
	int t_counts;
	int g_counts;
	int c_counts;
	int n_counts;





	//gzip random access index
	zran_index_t *gzip_index;

} FastaObject;

extern PyTypeObject pyfastx_FastaType;

PyObject* build_index(FastxObject *self, PyObject *args, PyObject *kwargs);
PyObject* test(FastxObject *self, PyObject *args);
PyObject* fastx_tp_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);
PyObject* fastx_tp_iter(FastxObject *self);
PyObject* fastx_tp_next(FastxObject *self);
int fastx_get_item(FastxObject *self, PyObject *key);
PyObject *fastx_get_key(FastxObject *self, PyObject *key);
void fastx_tp_dealloc(FastxObject *self);
int fastx_get_len(FastxObject *self);
int fastx_get_val(FastxObject *self, PyObject *key, PyObject *val);

#endif