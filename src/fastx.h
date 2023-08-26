#ifndef PYFASTX_FASTX_H
#define PYFASTX_FASTX_H
#include <Python.h>
#include "zlib.h"
#include "kseq.h"

typedef struct {
	PyObject_HEAD

	//fasta or fastq file path and name
	PyObject* file_obj;

	//always output uppercase sequence
	int uppercase;

	//file format fasta or fastq
	int format;

	//output comment or not
	int comment;

	//gzip open file handle
	gzFile gzfd;

	//kseqs for reading from fasta/q
	kseq_t* kseqs;

	PyObject* (*func) (kseq_t *);

} pyfastx_Fastx;

extern PyTypeObject pyfastx_FastxType;

#endif