#ifndef PYFASTX_UTIL_H
#define PYFASTX_UTIL_H
#include <Python.h>
int file_exists(char *file_name);
void upper_string(char *str);
int is_gzip(FILE *fd);

PyObject* clean_seq(PyObject *self, PyObject *args);
PyObject* sub_seq(PyObject *self, PyObject *args);

#endif