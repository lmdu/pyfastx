#ifndef PYFASTX_UTIL_H
#define PYFASTX_UTIL_H
#include <Python.h>
#include <stdint.h>
#include "sqlite3.h"
#include "zran.h"

uint16_t file_exists(char *file_name);
void remove_space(char *str);
void upper_string(char *str);
void remove_space_uppercase(char *str);
void reverse_seq(char *seq);
void reverse_complement_seq(char *seq);

int is_gzip_format(char *file_name);
//void truncate_seq(char *seq, uint32_t start, uint32_t end);
void complement_seq(char *seq);
void reverse_seq(char *seq);
uint32_t sum_array(uint32_t arr[], int num);
//char *int_to_str(int c);
int is_subset(char *seq1, char *seq2);
//PyObject* make_large_sequence(char *seq);
//int integer_check(PyObject* num);
//int64_t integer_to_long(PyObject* num);

//PyObject* clean_seq(PyObject *self, PyObject *args);
//PyObject* sub_seq(PyObject *self, PyObject *args);

//int64_t zran_readline(zran_index_t *index, char *linebuf, uint32_t bufsize);
void pyfastx_build_gzip_index(zran_index_t* gzip_index, sqlite3* index_db, char* index_file);
void pyfastx_load_gzip_index(zran_index_t* gzip_index, sqlite3* index_db, char* index_file);

#endif