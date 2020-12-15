#define PY_SSIZE_T_CLEAN
#include "fastx.h"
#include "util.h"

PyObject *pyfastx_fastx_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	//fasta or fastq file path
	Py_ssize_t file_len;
	char *file_name;
	int uppercase = 0;
	char *format = "auto";

	pyfastx_Fastx *obj;

	static char* keywords[] = {"file_name", "format", "uppercase", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|si", keywords, &file_name, &file_len, &format, &uppercase)) {
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

	//initial kseq
	obj->kseqs = kseq_init(obj->gzfd);

	return (PyObject *)obj;
}

void pyfastx_fastx_dealloc(pyfastx_Fastx *self) {
	free(self->file_name);
	kseq_destroy(self->kseqs);
	gzclose(self->gzfd);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *pyfastx_fastx_fasta(kseq_t* kseqs) {
	return Py_BuildValue("sss", kseqs->name.s, kseqs->seq.s, kseqs->comment.s);
}

PyObject *pyfastx_fastx_fasta_upper(kseq_t* kseqs) {
	upper_string(kseqs->seq.s, kseqs->seq.l);
	return pyfastx_fastx_fasta(kseqs);
}

PyObject *pyfastx_fastx_fastq(kseq_t* kseqs) {
	return Py_BuildValue("ssss", kseqs->name.s, kseqs->seq.s, kseqs->qual.s, kseqs->comment.s);
}

PyObject *pyfastx_fastx_iter(pyfastx_Fastx *self) {
	gzrewind(self->gzfd);

	if (self->format == 1) {
		if (self->uppercase) {
			self->func = pyfastx_fastx_fasta_upper;
		} else {
			self->func = pyfastx_fastx_fasta;
		}
	} else {
		self->func = pyfastx_fastx_fastq;
	}

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
	return PyUnicode_FromFormat("<Fastx> iterator for %s", self->file_name);
}

PyTypeObject pyfastx_FastxType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Fastx",                        /* tp_name */
    sizeof(pyfastx_Fastx),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_fastx_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)pyfastx_fastx_repr,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                   /* tp_as_sequence */
    0,                   /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,             /* tp_flags */
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    (getiterfunc)pyfastx_fastx_iter,     /* tp_iter */
    (iternextfunc)pyfastx_fastx_next,    /* tp_iternext */
    0,          /* tp_methods */
    0,          /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    pyfastx_fastx_new,              /* tp_new */
};