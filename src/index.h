#ifndef PYFASTX_INDEX_H
#define PYFASTX_INDEX_H
#include "Python.h"
#include "sqlite3.h"
#include "zlib.h"
#include "zran.h"
#include "kseq.h"

typedef struct {
	PyObject_HEAD

	//index file path
	char* index_file;

	//always output uppercase
	int uppercase;

	//is gzip compressed file
	//0 not gzip file
	//1 is gzip file
	int gzip_format;

	//open file handle
	FILE* fd;

	//gzip open file handle
	gzFile gzfd;
	
	//kseqs for reading from fasta
	kseq_t* kseqs;

	//sqlite3 handle for index
	sqlite3* index_db;

	//gzip random access index
	zran_index_t* gzip_index;

	//cache name
	char* cache_name;

	//cache start and end position
	int cache_start;
	int cache_end;

	//cache sequence
	char* cache_seq;

} pyfastx_Index;


void pyfastx_build_gzip_index(pyfastx_Index *self);
void pyfastx_load_gzip_index(pyfastx_Index *self);
void pyfastx_create_index(pyfastx_Index *self);
void pyfastx_load_index(pyfastx_Index *self);
void pyfastx_build_index(pyfastx_Index *self);
void pyfastx_rewind_index(pyfastx_Index *index);
void pyfastx_index_free(pyfastx_Index *self);

PyObject *pyfastx_get_next_seq(pyfastx_Index *index);
PyObject *pyfastx_index_make_seq(pyfastx_Index *self, sqlite3_stmt *stmt);
PyObject *pyfastx_index_get_seq_by_name(pyfastx_Index *self, char *name);
PyObject *pyfastx_index_get_seq_by_id(pyfastx_Index *self, int id);

pyfastx_Index *pyfastx_init_index(char* file_path, int uppercase);
char *pyfastx_index_get_sub_seq(pyfastx_Index *self, char *name, int64_t offset, int bytes, int start, int end, int normal);
char *pyfastx_index_get_full_seq(pyfastx_Index *self, char *name);

#endif