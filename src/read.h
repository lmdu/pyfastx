#ifndef PYFASTX_READ_H
#define PYFASTX_READ_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "zran.h"
#include "kseq.h"
#include "fastq.h"
#include "util.h"
#include "sqlite3.h"

typedef struct {
	PyObject_HEAD

	//read order id
	Py_ssize_t id;

	//read length
	Py_ssize_t read_len;

	//description length
	int desc_len;

	//seq start offset
	Py_ssize_t seq_offset;

	//quality start offset
	Py_ssize_t qual_offset;

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