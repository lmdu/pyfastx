#include "read.h"
#include "util.h"
#include "time.h"
#include "structmember.h"

void pyfastx_read_dealloc(pyfastx_Read *self) {
    free(self->name);

    if (self->seq) {
        free(self->seq);
    }

    if (self->qual) {
        free(self->qual);
    }

    if (self->raw) {
        free(self->raw);
    }

    if (self->desc) {
        free(self->desc);
    }

    Py_TYPE(self)->tp_free((PyObject *)self);
}

int pyfastx_read_length(pyfastx_Read *self) {
    return self->read_len;
}

//read content from buff or file
void pyfastx_read_continue_reader(pyfastx_Read *self) {
    int32_t slice_offset;
    int32_t slice_length;
    int32_t residue_len;
    int32_t read_len;
    int32_t cache_len;
    int64_t offset;

    //read raw string offset
    offset = self->seq_offset - self->desc_len - 1;

    //read raw string length
    residue_len = self->qual_offset + self->read_len - offset + 1;
    read_len = 0;

    self->raw = (char *)malloc(residue_len + 2);

    while (residue_len > 0) {
        if (offset >= self->fastq->cache_soff && offset < self->fastq->cache_eoff) {
            slice_offset = offset - self->fastq->cache_soff;
            slice_length = self->fastq->cache_eoff - offset;

            if (slice_length >= residue_len) {
                cache_len = residue_len;
            } else {
                cache_len = self->fastq->cache_eoff - self->fastq->cache_soff;
            }
            
            memcpy(self->raw+read_len, self->fastq->cache_buff+slice_offset, cache_len);
            read_len += cache_len;
            residue_len -= cache_len;
        }

        if (residue_len > 0) {
            self->fastq->cache_soff = self->fastq->cache_eoff;
            gzread(self->fastq->gzfd, self->fastq->cache_buff, 1048576);
            self->fastq->cache_eoff = gztell(self->fastq->gzfd);
        }
    }

    self->desc = (char *)malloc(self->desc_len + 1);
    memcpy(self->desc, self->raw, self->desc_len);

    if (self->raw[read_len-1] == '\r') {
        self->raw[read_len] = '\n';
        self->raw[read_len+1] = '\0';
        self->desc[self->desc_len-1] = '\0';
    } else {
        self->raw[read_len] = '\0';
        self->desc[self->desc_len] = '\0';
    }

    self->seq = (char *)malloc(self->read_len + 1);
    memcpy(self->seq, self->raw + self->seq_offset - offset, self->read_len);
    self->seq[self->read_len] = '\0';

    self->qual = (char *)malloc(self->read_len + 1);
    memcpy(self->qual, self->raw + self->qual_offset - offset, self->read_len);
    self->qual[self->read_len] = '\0';
}

void pyfastx_read_random_reader(pyfastx_Read *self, char *buff, int64_t offset, uint32_t bytes) {
    if (self->fastq->gzip_format) {
        zran_seek(self->fastq->gzip_index, offset, SEEK_SET, NULL);
        zran_read(self->fastq->gzip_index, buff, bytes);
    } else {
        gzseek(self->fastq->gzfd, offset, SEEK_SET);
        gzread(self->fastq->gzfd, buff, bytes);
    }
}

PyObject* pyfastx_read_raw(pyfastx_Read *self, void* closure) {
    if (self->fastq->iterating) {
        pyfastx_read_continue_reader(self);
    }
    else if (!self->raw) {
        int64_t new_offset;
        int64_t new_bytelen;
        new_offset = self->seq_offset - self->desc_len - 1;
        new_bytelen = self->qual_offset + self->read_len - new_offset + 1;

        self->raw = (char *)malloc(new_bytelen + 2);

        pyfastx_read_random_reader(self, self->raw, new_offset, new_bytelen);

        if (self->raw[new_bytelen-1] == '\r') {
            self->raw[new_bytelen] = '\n';
            self->raw[new_bytelen+1] = '\0';
        } else {
            self->raw[new_bytelen] = '\0';
        }
    }

    return Py_BuildValue("s", self->raw);
}

PyObject* pyfastx_read_seq(pyfastx_Read *self, void* closure) {
    if (self->fastq->iterating) {
        pyfastx_read_continue_reader(self);
    }
    else if (!self->seq) {
        self->seq = (char *)malloc(self->read_len + 1);
        pyfastx_read_random_reader(self, self->seq, self->seq_offset, self->read_len);
        self->seq[self->read_len] = '\0';
    }

    return Py_BuildValue("s", self->seq);
}

PyObject* pyfastx_read_description(pyfastx_Read *self, void* closure) {
    if (self->fastq->iterating) {
        pyfastx_read_continue_reader(self);
    }
    else if (!self->desc) {
        int64_t new_offset;

        new_offset = self->seq_offset - self->desc_len - 1;
        self->desc = (char *)malloc(self->desc_len + 1);

        pyfastx_read_random_reader(self, self->desc, new_offset, self->desc_len);

        if (self->desc[self->desc_len-1] == '\r') {
            self->desc[self->desc_len-1] = '\0';
        } else {
            self->desc[self->desc_len] = '\0';
        }
    }

    return Py_BuildValue("s", self->desc);
}

PyObject* pyfastx_read_qual(pyfastx_Read *self, void* closure) {
    if (self->fastq->iterating) {
        pyfastx_read_continue_reader(self);
    }
    else if (!self->qual) {
        self->qual = (char *)malloc(self->read_len + 1);
        pyfastx_read_random_reader(self, self->qual, self->qual_offset, self->read_len);
        self->qual[self->read_len] = '\0';
    }

    return Py_BuildValue("s", self->qual);
}

PyObject* pyfastx_read_quali(pyfastx_Read *self, void* closure) {
    int phred;
    int i;
    PyObject *quals;
    PyObject *q;

    if (!self->qual) {
        pyfastx_read_qual(self, NULL);
    }

    phred = self->fastq->phred ? self->fastq->phred : 33;

    quals = PyList_New(0);
    for (i = 0; i < self->read_len; i++) {
        q = Py_BuildValue("i", self->qual[i] - phred);
        PyList_Append(quals, q);
        Py_DECREF(q);
    }

    return quals;
}

PyObject* pyfastx_read_repr(pyfastx_Read *self) {
    return PyUnicode_FromFormat("<Read> %s with length of %d", self->name, self->read_len);
}

PyObject* pyfastx_read_str(pyfastx_Read *self) {
    return pyfastx_read_seq(self, NULL);
}

static PyGetSetDef pyfastx_read_getsets[] = {
    {"raw", (getter)pyfastx_read_raw, NULL, NULL, NULL},
    {"seq", (getter)pyfastx_read_seq, NULL, NULL, NULL},
    {"qual", (getter)pyfastx_read_qual, NULL, NULL, NULL},
    {"quali", (getter)pyfastx_read_quali, NULL, NULL, NULL},
    {"description", (getter)pyfastx_read_description, NULL, NULL, NULL},
    {NULL}
};

static PyMappingMethods pyfastx_read_as_mapping = {
	(lenfunc)pyfastx_read_length,
	//(binaryfunc)pyfastx_fasta_subscript,
	0,
};

static PyMemberDef pyfastx_read_members[] = {
    {"id", T_ULONGLONG, offsetof(pyfastx_Read, id), READONLY},
	{"name", T_STRING, offsetof(pyfastx_Read, name), READONLY},
	//{"size", T_LONG, offsetof(pyfastx_Read, seq_length), READONLY},
	//{"gc_content", T_FLOAT, offsetof(pyfastx_Read, gc_content), READONLY},
	//{"composition", T_OBJECT, offsetof(pyfastx_Read, composition), READONLY},
	{NULL}
};

PyTypeObject pyfastx_ReadType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Read",                        /* tp_name */
    sizeof(pyfastx_Read),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_read_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)pyfastx_read_repr,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                   /* tp_as_sequence */
    &pyfastx_read_as_mapping,                   /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    (reprfunc)pyfastx_read_str,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,             /* tp_flags */
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,     /* tp_iter */
    0,    /* tp_iternext */
    0,          /* tp_methods */
    pyfastx_read_members,          /* tp_members */
    pyfastx_read_getsets,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    PyType_GenericNew,              /* tp_new */
};