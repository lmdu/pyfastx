#ifndef PYFASTX_INDEX_H
#define PYFASTX_INDEX_H
#include "Python.h"
#include <stdint.h>
#include "sqlite3.h"
#include "zlib.h"
#include "zran.h"
#include "kseq.h"

typedef struct {
	PyObject_HEAD

	//index file path
	char* index_file;

	//always output uppercase
	uint8_t uppercase;

	//full name
	uint8_t full_name;

	//is gzip compressed file
	//0 not gzip file
	//1 is gzip file
	uint8_t gzip_format;

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

	//cahce seq id
	uint32_t cache_chrom;

	//cache seq start and end position
	uint32_t cache_start;
	uint32_t cache_end;

	//cache seq is complete or not
	uint8_t cache_full;

	//cache seq name string
	kstring_t cache_name;

	//cache real seq
	kstring_t cache_seq;

	//key function
	PyObject* key_func;

	//enter iterating loop
	uint8_t iterating;

	//prepared sql
	sqlite3_stmt *iter_stmt;
	sqlite3_stmt *uid_stmt;
	sqlite3_stmt *seq_stmt;

} pyfastx_Index;

//void pyfastx_build_gzip_index(pyfastx_Index *self);
//void pyfastx_load_gzip_index(pyfastx_Index *self);
void pyfastx_create_index(pyfastx_Index *self);
void pyfastx_load_index(pyfastx_Index *self);
void pyfastx_build_index(pyfastx_Index *self);
void pyfastx_rewind_index(pyfastx_Index *index);
void pyfastx_index_free(pyfastx_Index *self);
void pyfastx_index_cache_clear(pyfastx_Index *self);
//void pyfastx_index_continue_read(pyfastx_Index *self, char *buff, int64_t offset, uint32_t bytes);

PyObject *pyfastx_index_next_seq(pyfastx_Index *self);
PyObject *pyfastx_index_next_upper_seq(pyfastx_Index *self);
PyObject *pyfastx_index_next_full_name_seq(pyfastx_Index *self);
PyObject *pyfastx_index_next_full_name_upper_seq(pyfastx_Index *self);
PyObject *pyfastx_index_next_with_index_seq(pyfastx_Index *self);

PyObject *pyfastx_index_make_seq(pyfastx_Index *self, sqlite3_stmt *stmt);
PyObject *pyfastx_index_get_seq_by_name(pyfastx_Index *self, PyObject *name);
PyObject *pyfastx_index_get_seq_by_id(pyfastx_Index *self, uint32_t id);

pyfastx_Index *pyfastx_init_index(char* file_path, int file_len, int uppercase, int full_name, int memory_index, PyObject* key_func);
//char *pyfastx_index_get_sub_seq(pyfastx_Index *self, pyfastx_Sequence *seq);
//char *pyfastx_index_get_full_seq(pyfastx_Index *self, uint32_t chrom);
void pyfastx_index_random_read(pyfastx_Index* self, char* buff, int64_t offset, uint32_t bytes);
void pyfastx_index_fill_cache(pyfastx_Index* self, int64_t offset, uint32_t size);

#endif