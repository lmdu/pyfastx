/*
 * zran_file_util.c - file utilities used in zran.c to manipulate either
 * Python file-like objects or file descriptors.
 *
 * Functions which interact with Python file-likes will acquire and release
 * the GIL as needed.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "zran_file_util.h"

#ifdef _WIN32
#define FSEEK _fseeki64
#define FTELL _ftelli64
#include "windows.h"
#include "io.h"
#else
#define FSEEK fseeko
#define FTELL ftello
#endif

/*
 * The zran functions are typically called with the GIL released. These
 * macros are used to temporarily (re-)acquire and release the GIL when
 * interacting with Python file-like objects.
 */
#define _ZRAN_FILE_UTIL_ACQUIRE_GIL \
    PyGILState_STATE s;             \
    s = PyGILState_Ensure();

#define _ZRAN_FILE_UTIL_RELEASE_GIL \
    PyGILState_Release(s);


/*
 * Implements a method analogous to fread that is performed on Python
 * file-like objects.
 */
size_t _fread_python(void *ptr, size_t size, size_t nmemb, PyObject *f) {

    PyObject  *data = NULL;
    char      *buf;
    Py_ssize_t len;

    _ZRAN_FILE_UTIL_ACQUIRE_GIL

    if ((data = PyObject_CallMethod(f, "read", "(n)", size * nmemb)) == NULL)
        goto fail;
    if ((buf = PyBytes_AsString(data)) == NULL)
        goto fail;
    if ((len = PyBytes_Size(data)) == -1)
        goto fail;

    memmove(ptr, buf, (size_t) len);

    Py_DECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return (size_t) len / size;

fail:
    Py_XDECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return 0;
}

/*
 * Implements a method analogous to ftell that is performed on Python
 * file-like objects.
 */
int64_t _ftell_python(PyObject *f) {
    PyObject *data = NULL;
    int64_t   result;

    _ZRAN_FILE_UTIL_ACQUIRE_GIL

    data = PyObject_CallMethod(f, "tell", NULL);
    if (data == NULL)
        goto fail;

    result = PyLong_AsLong(data);
    if (result == -1 && PyErr_Occurred())
        goto fail;

    Py_DECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return result;

fail:
    Py_XDECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return -1;
}

/*
 * Implements a method analogous to fseek that is performed on Python
 * file-like objects.
 */
int _fseek_python(PyObject *f, int64_t offset, int whence) {

    PyObject *data = NULL;

    _ZRAN_FILE_UTIL_ACQUIRE_GIL

    /*
     * The seek method returns type long, which has
     * different sizes on different platforms
     */
    if (sizeof(long) == 8)
        data = PyObject_CallMethod(f, "seek", "(l,i)", offset, whence);
    else if (sizeof(long long) == 8)
        data = PyObject_CallMethod(f, "seek", "(L,i)", offset, whence);
    else
        goto fail;
    if (data == NULL)
        goto fail;

    Py_DECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return 0;

fail:
    Py_XDECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return -1;
}

/*
 * Implements a method analogous to feof that is performed on Python file-like
 * objects.  If f_ret, the number of bytes returned by the last read, is zero,
 * then we're at EOF.
 */
int _feof_python(PyObject *f, size_t f_ret) {
    return f_ret == 0;
}

/*
 * Implements a method analogous to ferror that is performed on Python
 * file-like objects.
 */
int _ferror_python(PyObject *f) {
    PyObject *result;

    _ZRAN_FILE_UTIL_ACQUIRE_GIL
    result = PyErr_Occurred();
    _ZRAN_FILE_UTIL_RELEASE_GIL

    if (result != NULL) return 1;
    else                return 0;
}

/*
 * Implements a method analogous to fflush that is performed on Python
 * file-like objects.
 */
int _fflush_python(PyObject *f) {
    PyObject *data = NULL;

    _ZRAN_FILE_UTIL_ACQUIRE_GIL
    if ((data = PyObject_CallMethod(f, "flush", NULL)) == NULL) goto fail;

    Py_DECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return 0;

fail:
    Py_XDECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return -1;
}

/*
 * Implements a method analogous to fwrite that is performed on Python
 * file-like objects.
 */
size_t _fwrite_python(const void *ptr,
                      size_t      size,
                      size_t      nmemb,
                      PyObject   *f) {

    PyObject *input = NULL;
    PyObject *data = NULL;
    long      len;

    _ZRAN_FILE_UTIL_ACQUIRE_GIL
    if ((input = PyBytes_FromStringAndSize(ptr, size * nmemb)) == NULL)
        goto fail;
    if ((data = PyObject_CallMethod(f, "write", "(O)", input)) == NULL)
        goto fail;

    #if PY_MAJOR_VERSION >= 3
    if ((len = PyLong_AsLong(data)) == -1 && PyErr_Occurred())
        goto fail;
    #else
    // In Python 2, a file object's write() method does not return the number
    // of bytes written, so let's just assume that everything has been written
    // properly.
    len = size * nmemb;
    #endif

    Py_DECREF(input);
    Py_DECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return (size_t) len / size;

fail:
    Py_XDECREF(input);
    Py_XDECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return 0;
}

/*
 * Implements a method analogous to getc that is performed on Python file-like
 * objects.
 */
int _getc_python(PyObject *f) {
    char buf;
    if (_fread_python(&buf, 1, 1, f) == 0) {
        // Reached EOF, or an error (in which case the error indicator is set).
        // Either way, we should return -1.
        return -1;
    }
    return buf;
}

/*
 * Calls the .seekable() method on Python file-like objects.
 */
int _seekable_python(PyObject *f) {
    PyObject *data = NULL;
    int64_t   result;

    _ZRAN_FILE_UTIL_ACQUIRE_GIL

    data = PyObject_CallMethod(f, "seekable", NULL);
    if (data == NULL)
        goto fail;

    result = PyLong_AsLong(data);
    if (result == -1 && PyErr_Occurred())
        goto fail;

    Py_DECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return result;

fail:
    Py_XDECREF(data);
    _ZRAN_FILE_UTIL_RELEASE_GIL
    return -1;
}


/*
 * Hacky implementation of seekable() for Python 2.
 * If f.tell() returns an error, assume that the
 * object is unseekable.
 */
int _seekable_python2(PyObject *f) {
    int64_t ret;
    ret = _ftell_python(f);
    if (ret < 0) {
        PyErr_Clear();
    }
    return ret >= 0;
}

/*
 * Calls ferror on fd if specified, otherwise the Python-specific method on f.
 */
int ferror_(FILE *fd, PyObject *f) {
    return fd != NULL ? ferror(fd) : _ferror_python(f);
}

/*
 * Calls fseek on fd if specified, otherwise the Python-specific method on f.
 */
int fseek_(FILE *fd, PyObject *f, int64_t offset, int whence) {
    return fd != NULL
        ? FSEEK(fd, offset, whence)
        : _fseek_python(f, offset, whence);
}

/*
 * Calls ftell on fd if specified, otherwise the Python-specific method on f.
 */
int64_t ftell_(FILE *fd, PyObject *f) {
    return fd != NULL ? FTELL(fd) : _ftell_python(f);
}

/*
 * Calls fread on fd if specified, otherwise the Python-specific method on f.
 */
size_t fread_(void *ptr, size_t size, size_t nmemb, FILE *fd, PyObject *f) {

  return fd != NULL
      ? fread(ptr, size, nmemb, fd)
      : _fread_python(ptr, size, nmemb, f);
}

/*
 * Calls feof on fd if specified, otherwise the Python-specific method on f.
 * If fd is not specified, requires f_ret, the number of bytes read on the last
 * read, to determine if the file is at EOF.
 */
int feof_(FILE *fd, PyObject *f, size_t f_ret) {
    return fd != NULL ? feof(fd): _feof_python(f, f_ret);
}

/*
 * Calls fflush on fd if specified, otherwise the Python-specific method on f.
 */
int fflush_(FILE *fd, PyObject *f) {
    return fd != NULL ? fflush(fd): _fflush_python(f);
}

/*
 * Calls fwrite on fd if specified, otherwise the Python-specific method on f.
 */
size_t fwrite_(const void *ptr,
               size_t      size,
               size_t      nmemb,
               FILE       *fd,
               PyObject   *f) {
    return fd != NULL
        ? fwrite(ptr, size, nmemb, fd)
        : _fwrite_python(ptr, size, nmemb, f);
}

/*
 * Calls getc on fd if specified, otherwise the Python-specific method on f.
 */
int getc_(FILE *fd, PyObject *f) {
    return fd != NULL ? getc(fd): _getc_python(f);
}

/*
 * Returns whether the given file is seekable. If fd is specified, assumes it's always seekable.
 * If f is specified, calls f.seekable() to see if the Python file object is seekable.
 */
int seekable_(FILE *fd, PyObject *f) {
    #if PY_MAJOR_VERSION > 2
    return fd != NULL ? 1: _seekable_python(f);
    #else
    return fd != NULL ? 1: _seekable_python2(f);
    #endif
}
