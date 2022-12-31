#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "fastx.h"
#include "util.h"

PyObject *pyfastx_fastx_fasta(kseq_t* kseqs) {
	return Py_BuildValue("ss", kseqs->name.s, kseqs->seq.s);
}

PyObject *pyfastx_fastx_fasta_comment(kseq_t* kseqs) {
	return Py_BuildValue("sss#", kseqs->name.s, kseqs->seq.s, kseqs->comment.s, kseqs->comment.l);
}

PyObject *pyfastx_fastx_fasta_upper(kseq_t* kseqs) {
	return pyfastx_fastx_fasta(kseqs);
}

PyObject *pyfastx_fastx_fasta_upper_comment(kseq_t* kseqs) {
	upper_string(kseqs->seq.s, kseqs->seq.l);
	return pyfastx_fastx_fasta_comment(kseqs);
}

PyObject *pyfastx_fastx_fastq(kseq_t* kseqs) {
	return Py_BuildValue("sss", kseqs->name.s, kseqs->seq.s, kseqs->qual.s);
}

PyObject *pyfastx_fastx_fastq_comment(kseq_t* kseqs) {
	return Py_BuildValue("ssss#", kseqs->name.s, kseqs->seq.s, kseqs->qual.s, kseqs->comment.s, kseqs->comment.l);
}

PyObject *pyfastx_fastx_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	int uppercase = 0;
	int comment = 0;

	char *file_name;
	char *format = "auto";

	Py_ssize_t file_len;

	pyfastx_Fastx *obj;

	static char* keywords[] = {"file_name", "format", "uppercase", "comment", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|sii", keywords, &file_name, &file_len, &format, &uppercase, &comment)) {
		return NULL;
	}

	if (!file_exists(file_name)) {
		PyErr_Format(PyExc_FileExistsError, "the input file %s does not exists", file_name);
		return NULL;
	}

	obj = (pyfastx_Fastx *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	obj->file_name = (char *)malloc(file_len + 1);
	strcpy(obj->file_name, file_name);

	//open the sequence file
	obj->gzfd = gzopen(file_name, "rb");

	//set file format
	if (strcmp(format, "auto") == 0) {
		obj->format = fasta_or_fastq(obj->gzfd);
	} else if (strcmp(format, "fasta") == 0) {
		obj->format = 1;
	} else if (strcmp(format, "fastq") == 0) {
		obj->format = 2;
	} else {
		obj->format = 0;
	}

	if (obj->format == 0) {
		PyErr_Format(PyExc_RuntimeError, "%s is not fasta or fastq sequence file", file_name);
		return NULL;
	}

	obj->uppercase = uppercase;
	obj->comment = comment;

	//initial kseq
	gzrewind(obj->gzfd);
	obj->kseqs = kseq_init(obj->gzfd);

	if (obj->format == 1) {
		if (obj->uppercase) {
			if (obj->comment) {
				obj->func = pyfastx_fastx_fasta_upper_comment;
			} else {
				obj->func = pyfastx_fastx_fasta_upper;
			}
		} else {
			if (obj->comment) {
				obj->func = pyfastx_fastx_fasta_comment;
			} else {
				obj->func = pyfastx_fastx_fasta;
			}
		}
	} else {
		if (obj->comment) {
			obj->func = pyfastx_fastx_fastq_comment;
		} else {
			obj->func = pyfastx_fastx_fastq;
		}
	}

	return (PyObject *)obj;
}

void pyfastx_fastx_dealloc(pyfastx_Fastx *self) {
	free(self->file_name);
	kseq_destroy(self->kseqs);
	gzclose(self->gzfd);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *pyfastx_fastx_iter(pyfastx_Fastx *self) {
	gzrewind(self->gzfd);
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fastx_next(pyfastx_Fastx *self) {
	if (kseq_read(self->kseqs) >= 0) {
		return self->func(self->kseqs);
	}

	return NULL;
}

PyObject *pyfastx_fastx_repr(pyfastx_Fastx *self) {
	if (self->format == 1) {
		return PyUnicode_FromFormat("<Fastx> fasta %s", self->file_name);
	} else {
		return PyUnicode_FromFormat("<Fastx> fastq %s", self->file_name);
	}
}

PyTypeObject pyfastx_FastxType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Fastx",
    .tp_basicsize = sizeof(pyfastx_Fastx),
    .tp_dealloc = (destructor)pyfastx_fastx_dealloc,
    .tp_repr = (reprfunc)pyfastx_fastx_repr,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_iter = (getiterfunc)pyfastx_fastx_iter,
    .tp_iternext = (iternextfunc)pyfastx_fastx_next,
    .tp_new = pyfastx_fastx_new,
};
