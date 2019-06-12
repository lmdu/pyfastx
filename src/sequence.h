#ifndef PYFASTX_SEQ_H
#define PYFASTX_SEQ_H
#include "fasta.h"

//make sequence class
typedef struct {
	PyObject_HEAD
	
	//sequence name
	char* name;

	//start position
	int start;

	//end position
	int end;

	//sequence length
	int seq_len;

	//GC content
	float gc_content;

	//a dict for storing ATGCN counts
	PyObject* composition;

	pyfastx_Fasta* fasta;

} pyfastx_Sequence;

extern PyTypeObject pyfastx_SeqType;

#endif