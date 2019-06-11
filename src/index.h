#ifndef PYFASTX_INDEX_H
#define PYFASTX_INDEX_H
#include "kseq.h"
#include "fasta.h"

KSEQ_DECLARE(gzFile)

static void _pyfastx_build_gzip_index(pyfastx_Fasta *self);
static void _pyfastx_load_gzip_index(pyfastx_Fasta *self);
PyObject *pyfastx_build_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs);
PyObject *get_sub_seq(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs);

#endif