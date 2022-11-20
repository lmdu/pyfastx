#ifndef __ZRAN_FILE_UTIL_H__
#define __ZRAN_FILE_UTIL_H__

/*
 * File utilities used to manipulate either
 * Python file-like objects or file descriptors.
 */

#include <stdlib.h>
#include <stdint.h>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

/*
 * Implements a method analogous to fread that is performed on Python
 * file-like objects.
 */
size_t _fread_python(void *ptr, size_t size, size_t nmemb, PyObject *f);

/*
 * Implements a method analogous to ftell that is performed on Python
 * file-like objects.
 */
int64_t _ftell_python(PyObject *f);

/*
 * Implements a method analogous to fseek that is performed on Python
 * file-like objects.
 */
int _fseek_python(PyObject *f, int64_t offset, int whence);

/*
 * Implements a method analogous to feof that is performed on Python file-like
 * objects.
 */
int _feof_python(PyObject *f, size_t f_ret);

/*
 * Implements a method analogous to ferror that is performed on Python
 * file-like objects.
 */
int _ferror_python(PyObject *f);

/*
 * Implements a method analogous to fflush that is performed on Python
 * file-like objects.
 */
int _fflush_python(PyObject *f);

/*
 * Implements a method analogous to fwrite that is performed on Python
 * file-like objects.
 */
size_t _fwrite_python(const void *ptr, size_t size, size_t nmemb, PyObject *f);

/*
 * Implements a method analogous to getc that is performed on Python
 * file-like objects.
 */
int _getc_python(PyObject *f);

/*
 * Calls the .seekable() method on Python file-like objects.
 */
int _seekable_python(PyObject *f);

/*
 * Calls ferror on fd if specified, otherwise the Python-specific method on f.
 */
int ferror_(FILE *fd, PyObject *f);

/*
 * Calls fseek on fd if specified, otherwise the Python-specific method on f.
 */
int fseek_(FILE *fd, PyObject *f, int64_t offset, int whence);

/*
 * Calls ftell on fd if specified, otherwise the Python-specific method on f.
 */
int64_t ftell_(FILE *fd, PyObject *f);

/*
 * Calls fread on fd if specified, otherwise the Python-specific method on f.
 */
size_t fread_(void *ptr, size_t size, size_t nmemb, FILE *fd, PyObject *f);

/*
 * Calls feof on fd if specified, otherwise the Python-specific method on f.
 */
int feof_(FILE *fd, PyObject *f, size_t f_ret);

/*
 * Calls fflush on fd if specified, otherwise the Python-specific method on f.
 */
int fflush_(FILE *fd, PyObject *f);

/*
 * Calls fwrite on fd if specified, otherwise the Python-specific method on f.
 */
size_t fwrite_(
    const void *ptr, size_t size, size_t nmemb, FILE *fd, PyObject *f);

/*
 * Calls getc on fd if specified, otherwise the Python-specific method on f.
 */
int getc_(FILE *fd, PyObject *f);

/*
 * Returns whether the given file is seekable. If fd is specified, assumes it's always seekable.
 * If f is specified, calls f.seekable() to see if the Python file object is seekable.
 */
int seekable_(FILE *fd, PyObject *f);

#endif /* __ZRAN_FILE_UTIL_H__ */
