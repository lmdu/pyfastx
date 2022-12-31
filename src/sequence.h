#ifndef PYFASTX_SEQ_H
#define PYFASTX_SEQ_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "index.h"
#include "kseq.h"
#include "util.h"

//make sequence class
typedef struct {
	PyObject_HEAD

	//sequence order
	Py_ssize_t id;
	
	//sequence name
	char* name;

	//sequence description
	char* desc;

	//raw string
	char* raw;

	//iter line
	kstring_t line;

	//start position
	Py_ssize_t start;

	//end position
	Py_ssize_t end;

	//sequence real length
	Py_ssize_t seq_len;

	//description length
	int desc_len;

	//subsequence parent length
	//uint32_t parent_len;

	//fasta index
	pyfastx_Index* index;

	//start offset in fasta file
	Py_ssize_t offset;

	//byte length for sequence
	Py_ssize_t byte_len;

	//each line length of sequence
	Py_ssize_t line_len;

	//line end length, \n or \r\n
	int end_len;

	//standard fasta format with same line length
	int normal;

	//complete sequence or subsequence
	int complete;

	//line iteration cache
	char* line_cache;

	//current cache position
	char* cache_pos;

} pyfastx_Sequence;

extern PyTypeObject pyfastx_SequenceType;

Py_ssize_t pyfastx_sequence_length(pyfastx_Sequence* self);
char *pyfastx_sequence_acquire(pyfastx_Sequence* self);
int pyfastx_sequence_contains(pyfastx_Sequence *self, PyObject *key);

char *pyfastx_sequence_get_subseq(pyfastx_Sequence* self);
char *pyfastx_sequence_get_fullseq(pyfastx_Sequence* self);

PyObject *pyfastx_sequence_seq(pyfastx_Sequence* self, void* closure);
PyObject *pyfastx_sequence_reverse(pyfastx_Sequence* self, void* closure);
PyObject *pyfastx_sequence_complement(pyfastx_Sequence* self, void* closure);
PyObject *pyfastx_sequence_antisense(pyfastx_Sequence* self, void* closure);

PyObject *pyfastx_sequence_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);
PyObject *pyfastx_seqeunce_subscript(pyfastx_Sequence* self, PyObject* item);
PyObject *pyfastx_sequence_str(pyfastx_Sequence* self);
PyObject *pyfastx_sequence_repr(pyfastx_Sequence* self);
PyObject *pyfastx_sequence_iter(pyfastx_Sequence* self);
PyObject *pyfastx_sequnece_next(pyfastx_Sequence* self);

PyObject *pyfastx_sequence_search(pyfastx_Sequence *self, PyObject *args, PyObject *kwargs);

#endif