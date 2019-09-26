#ifndef PYFASTX_UTIL_H
#define PYFASTX_UTIL_H
#include <Python.h>
#include <stdint.h>

uint16_t file_exists(char *file_name);
void remove_space(char *str);
void upper_string(char *str);
uint16_t is_gzip_format(char *file_name);
void truncate_seq(char *seq, uint32_t start, uint32_t end);
void complement_seq(char *seq);
void reverse_seq(char *seq);

PyObject* clean_seq(PyObject *self, PyObject *args);
PyObject* sub_seq(PyObject *self, PyObject *args);

//int64_t zran_readline(zran_index_t *index, char *linebuf, uint32_t bufsize);

#endif