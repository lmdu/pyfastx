#ifndef PYFASTX_UTIL_H
#define PYFASTX_UTIL_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>
//#include "unistd.h"
#include "sqlite3.h"
#include "zran.h"
#include "zlib.h"
#include "time.h"

int file_exists(PyObject *file_obj);
void upper_string(char *str, Py_ssize_t len);
Py_ssize_t remove_space(char *str, Py_ssize_t len);
Py_ssize_t remove_space_uppercase(char *str, Py_ssize_t len);
void reverse_seq(char *seq);
void reverse_complement_seq(char *seq);

int is_gzip_format(PyObject *file_obj);
//void truncate_seq(char *seq, uint32_t start, uint32_t end);
void complement_seq(char *seq);
void reverse_seq(char *seq);
Py_ssize_t sum_array(Py_ssize_t arr[], int num);
//char *int_to_str(int c);
int is_subset(char *seq1, char *seq2);
//PyObject* make_large_sequence(char *seq);
//int integer_check(PyObject* num);
//int64_t integer_to_long(PyObject* num);
gzFile pyfastx_gzip_open(PyObject *path, const char *mode);

//PyObject* clean_seq(PyObject *self, PyObject *args);
//PyObject* sub_seq(PyObject *self, PyObject *args);
char *str_n_str(char *haystack, char *needle, Py_ssize_t len, Py_ssize_t size);

//int64_t zran_readline(zran_index_t *index, char *linebuf, uint32_t bufsize);
void pyfastx_build_gzip_index(zran_index_t* gzip_index, sqlite3* index_db);
void pyfastx_load_gzip_index(zran_index_t* gzip_index, sqlite3* index_db);

//a simple fasta/q validator
int fasta_validator(gzFile fd);
int fastq_validator(gzFile fd);
int fasta_or_fastq(gzFile fd);

//read line
/*ssize_t get_until_delim(char **buf, int delimiter, FILE *fp);
ssize_t get_line(char **buf, FILE *fp);*/

//sqlite3 compatable with python GIL
//extracted from apsw project
//https://github.com/rogerbinns/apsw/blob/07571365b6fbb25e2691071998526c351b04a04d/src/util.c
/* call where no error is returned */
#define PYFASTX_SQLITE_CALL(x) \
	do { Py_BEGIN_ALLOW_THREADS { x; } Py_END_ALLOW_THREADS ; } while(0)

#define test_time(x) \
	do { clock_t s, e; s=clock(); x; e=clock(); fprintf(stderr, "time: %.12f\n", (double)(e-s)/CLOCKS_PER_SEC); } while(0)

//support large fseek offset
#ifdef _WIN32
	#define FSEEK _fseeki64
	#define FTELL _ftelli64
#else
	#define FSEEK fseeko
	#define FTELL ftello
#endif

#ifndef _Py_CAST
#define _Py_CAST(type, expr) ((type)(expr))
#endif

#ifndef _PyObject_CAST
#define _PyObject_CAST(op) _Py_CAST(PyObject*, op)
#endif

#if PY_VERSION_HEX < 0x030A00A3 && !defined(Py_NewRef)
static inline PyObject* _Py_NewRef(PyObject *obj) {
	Py_INCREF(obj);
	return obj;
}
#define Py_NewRef(obj) _Py_NewRef(_PyObject_CAST(obj))
#endif

#if PY_VERSION_HEX < 0x030A00A3 && !defined(Py_XNewRef)
static inline PyObject* _Py_XNewRef(PyObject *obj) {
	Py_XINCREF(obj);
	return obj;
}
#define Py_XNewRef(obj) _Py_XNewRef(_PyObject_CAST(obj))
#endif

#endif