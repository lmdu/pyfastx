#ifndef PYFASTX_FASTA_H
#define PYFASTX_FASTA_H
#include "Python.h"
#include "zlib.h"
#include "kseq.h"
#include "index.h"

//make sequence iterator
typedef struct {
	PyObject_HEAD

	//fasta or fastq file path and name
	char* file_name;

	//always output uppercase sequence
	int uppercase;

	//total sequence counts
	int seq_counts;

	//total sequence length (bp)
	long seq_length;

	//a float for GC content (%)
	float gc_content;

	//a dict for storing A T G C N (unknown base) counts in fasta
	PyObject* composition;

	//index for fast random access to sequence
	pyfastx_Index* index;

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

void pyfastx_calc_fasta_attrs(pyfastx_Fasta *self);

PyObject *fasta_build_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs);
PyObject *fasta_rebuild_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs);
//PyObject *get_sub_seq(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs);

#endif