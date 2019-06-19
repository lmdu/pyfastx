#ifndef PYFASTX_SEQ_H
#define PYFASTX_SEQ_H
#include "Python.h"
#include "index.h"

//make sequence class
typedef struct {
	PyObject_HEAD
	
	//sequence name
	char* name;

	//start position
	int start;

	//end position
	int end;

	//sequence real length
	int seq_len;

	//GC content
	float gc_content;

	//a dict for storing ATGCN counts
	PyObject* composition;

	//fasta index
	pyfastx_Index* index;

	//id in db
	//int index_id;

	//start offset in fasta file
	int offset;

	//byte length for sequence
	int byte_len;

	//each line length of sequence
	int line_len;

	//line end length, \n or \r\n
	int end_len;

	//standard fasta format with same line length
	int normal;

} pyfastx_Sequence;

extern PyTypeObject pyfastx_SequenceType;

int pyfastx_sequence_length(pyfastx_Sequence* self);

PyObject *pyfastx_sequence_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);
PyObject *pyfastx_seqeunce_subscript(pyfastx_Sequence* self, PyObject* item);
PyObject *pyfastx_sequence_str(pyfastx_Sequence* self);

#endif