#ifndef AC_KSEQ_H
#define AC_KSEQ_H
#include "Python.h"
#include "zlib.h"

typedef struct __kstream_t {
	unsigned char *buf;
	int begin, end, is_eof;
	gzFile f;
} kstream_t;

typedef struct __kstring_t {
	size_t l, m;
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
int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append);
int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret);
kseq_t *kseq_init(gzFile fd);
void kseq_rewind(kseq_t *ks);
void kseq_destroy(kseq_t *ks);
int kseq_read(kseq_t *seq);

#endif