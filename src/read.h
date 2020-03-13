#ifndef PYFASTX_READ_H
#define PYFASTX_READ_H
#include "Python.h"
#include "zran.h"
#include "kseq.h"
#include "sqlite3.h"

typedef struct {
	PyObject_HEAD

	//read order id
	uint64_t id;

	//read name
	char* name;

	int read_len;

	int desc_len;

	//seq start offset
	int64_t seq_offset;

	//quality start offset
	int64_t qual_offset;

	//gzip index read handle
	gzFile gzfd;

	//file handle
	FILE* fd;

	//gzip index
	zran_index_t* gzip_index;

	//gzip format
	uint16_t gzip_format;

	//phred
	uint16_t phred;

	//seq and quality content
	char *seq;
	char *qual;

} pyfastx_Read;

extern PyTypeObject pyfastx_ReadType;

#endif