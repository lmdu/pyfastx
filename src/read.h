#ifndef PYFASTX_READ_H
#define PYFASTX_READ_H
#include "Python.h"
#include "zran.h"
#include "kseq.h"
#include "fastq.h"
#include "sqlite3.h"

typedef struct {
	PyObject_HEAD

	//read order id
	uint64_t id;

	//read length
	int read_len;

	//description length
	int desc_len;

	//seq start offset
	int64_t seq_offset;

	//quality start offset
	int64_t qual_offset;

	//parent fastq
	pyfastx_FastqMiddleware *middle;

	//read name
	char *name;

	//sequence
	char *seq;

	//quality
	char *qual;

	//raw string with name quality
	char *raw;

	//description
	char *desc;

} pyfastx_Read;

extern PyTypeObject pyfastx_ReadType;

#endif