#ifndef AC_KSEQ_H
#define AC_KSEQ_H
#define PY_SSIZE_T_CLEAN
#include <Python.h>
//#include <stdint.h>
#include <ctype.h>
#include "zlib.h"

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2
#define BUF_SIZE     1048576

#define ks_err(ks) ((ks)->end == -1)
#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

//kstring_t init
#define kstring_init(ks) ((ks).l = (ks).m = 0, (ks).s = NULL)

typedef struct __kstream_t {
	unsigned char *buf;
	//int64_t begin, end, is_eof;
	Py_ssize_t begin, end, is_eof;
	gzFile f;
} kstream_t;

typedef struct __kstring_t {
	//int64_t l, m;
	Py_ssize_t l, m;
	char *s;
} kstring_t;

typedef struct {
	kstring_t name, comment, seq, qual;
	int last_char;
	kstream_t *f;
} kseq_t;

kstream_t *ks_init(gzFile f);
void ks_destroy(kstream_t *ks);
int ks_getc(kstream_t *ks);
Py_ssize_t ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append);
Py_ssize_t ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret);
kseq_t *kseq_init(gzFile fd);
void kseq_rewind(kseq_t *ks);
void kseq_destroy(kseq_t *ks);
Py_ssize_t kseq_read(kseq_t *seq);

#endif