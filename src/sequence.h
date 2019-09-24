#ifndef PYFASTX_SEQ_H
#define PYFASTX_SEQ_H
#include "Python.h"
#include <stdint.h>
#include "index.h"

//make sequence class
typedef struct {
	PyObject_HEAD

	//sequence order
	uint32_t id;
	
	//sequence name
	char* name;

	//start position
	uint32_t start;

	//end position
	uint32_t end;

	//sequence real length
	uint32_t seq_len;

	//subsequence parent length
	uint32_t parent_len;

	//GC content
	float gc_content;

	//GC skew
	float gc_skew;

	//a dict for storing ATGCN counts
	PyObject* composition;

	//fasta index
	pyfastx_Index* index;

	//id in db
	//int index_id;

	//start offset in fasta file
	int64_t offset;

	//byte length for sequence
	uint32_t byte_len;

	//each line length of sequence
	uint32_t line_len;

	//line end length, \n or \r\n
	uint16_t end_len;

	//standard fasta format with same line length
	uint16_t normal;

	//for iteration
	kstream_t *ks;

} pyfastx_Sequence;

extern PyTypeObject pyfastx_SequenceType;

uint32_t pyfastx_sequence_length(pyfastx_Sequence* self);
char *pyfastx_sequence_acquire(pyfastx_Sequence* self);
uint16_t pyfastx_sequence_contains(pyfastx_Sequence *self, PyObject *key);

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