#define PY_SSIZE_T_CLEAN
#include <Python.h>
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

    Py_DECREF(self->middle->fastq);
    self->middle = NULL;

    Py_TYPE(self)->tp_free((PyObject *)self);
}

Py_ssize_t pyfastx_read_length(pyfastx_Read *self) {
    return self->read_len;
}

void pyfastx_read_random_reader(pyfastx_Read *self, char *buff, Py_ssize_t offset, Py_ssize_t bytes) {
    if (self->middle->gzip_format) {
        zran_seek(self->middle->gzip_index, offset, SEEK_SET, NULL);
        zran_read(self->middle->gzip_index, buff, bytes);
    } else {
        FSEEK(self->middle->fd, offset, SEEK_SET);
        fread(buff, bytes, 1, self->middle->fd);
    }
}

//read content from buff or file
void pyfastx_read_continue_reader(pyfastx_Read *self) {
    Py_ssize_t slice_offset;
    Py_ssize_t slice_length;
    Py_ssize_t residue_len;
    Py_ssize_t read_len;
    Py_ssize_t cache_len;
    Py_ssize_t offset;

    //read raw string offset
    offset = self->seq_offset - self->desc_len - 1;

    //read raw string length
    residue_len = self->qual_offset + self->read_len - offset + 2;
    read_len = 0;
    cache_len = 0;

    self->raw = (char *)malloc(residue_len + 1);

    if (offset < self->middle->cache_soff) {
        pyfastx_read_random_reader(self, self->raw, offset, residue_len);
    } else {
        while (residue_len > 0) {
            if (offset >= self->middle->cache_soff && offset < self->middle->cache_eoff) {
                slice_offset = offset - self->middle->cache_soff;
                slice_length = self->middle->cache_eoff - offset;

                cache_len = slice_length >= residue_len ? residue_len : slice_length;

                memcpy(self->raw+read_len, self->middle->cache_buff+slice_offset, cache_len);
                read_len += cache_len;
                residue_len -= cache_len;
            } else {
                self->middle->cache_soff = self->middle->cache_eoff;
                gzread(self->middle->gzfd, self->middle->cache_buff, 1048576);
                self->middle->cache_eoff = gztell(self->middle->gzfd);
            }
        }
    }

    self->desc = (char *)malloc(self->desc_len + 1);
    memcpy(self->desc, self->raw, self->desc_len);
    if (self->desc[self->desc_len-1] == '\r') {
        self->desc[self->desc_len-1] = '\0';
    } else {
       self->desc[self->desc_len] = '\0';
    }

    if (self->raw[read_len-2] == '\n') {
        self->raw[read_len-1] = '\0';
    } else if (self->raw[read_len-2] == '\r' && self->raw[read_len-1] == '\n') {
        self->raw[read_len] = '\0';
    } else {
        self->raw[read_len-2] = '\0';
    }

    self->seq = (char *)malloc(self->read_len + 1);
    memcpy(self->seq, self->raw + self->seq_offset - offset, self->read_len);
    self->seq[self->read_len] = '\0';

    self->qual = (char *)malloc(self->read_len + 1);
    memcpy(self->qual, self->raw + self->qual_offset - offset, self->read_len);
    self->qual[self->read_len] = '\0';
}

PyObject* pyfastx_read_raw(pyfastx_Read *self, void* closure) {
    Py_ssize_t new_offset;
    Py_ssize_t new_bytelen;

    if (! self->raw) {
        if (self->middle->iterating) {
            pyfastx_read_continue_reader(self);
        } else {
            new_offset = self->seq_offset - self->desc_len - 1;
            new_bytelen = self->qual_offset + self->read_len - new_offset + 2;

            self->raw = (char *)malloc(new_bytelen + 1);

            pyfastx_read_random_reader(self, self->raw, new_offset, new_bytelen);

            if (self->raw[new_bytelen-2] == '\n') {
                self->raw[new_bytelen-1] = '\0';
            } else if (self->raw[new_bytelen-2] == '\r' && self->raw[new_bytelen-1] == '\n') {
                self->raw[new_bytelen] = '\0';
            } else {
                self->raw[new_bytelen-2] = '\0';
            }
        }
    }

    return Py_BuildValue("s", self->raw);
}

void pyfastx_read_get_seq(pyfastx_Read* self) {
    if (! self->seq) {
        if (self->middle->iterating) {
            pyfastx_read_continue_reader(self);
        } else {
            self->seq = (char *)malloc(self->read_len + 1);
            pyfastx_read_random_reader(self, self->seq, self->seq_offset, self->read_len);
            self->seq[self->read_len] = '\0';
        }
    }
}

PyObject* pyfastx_read_seq(pyfastx_Read *self, void* closure) {
    pyfastx_read_get_seq(self);
    return Py_BuildValue("s", self->seq);
}

PyObject* pyfastx_read_reverse(pyfastx_Read *self, void* closure) {
    char *data;
    PyObject* ret;

    pyfastx_read_get_seq(self);

    ret = PyUnicode_New(self->read_len, 127);
    data = (char *)PyUnicode_1BYTE_DATA(ret);
    memcpy(data, self->seq, self->read_len);
 
    reverse_seq(data);

    return ret;
}

PyObject* pyfastx_read_complement(pyfastx_Read *self, void* closure) {
    char *data;
    PyObject* ret;

    pyfastx_read_get_seq(self);

    ret = PyUnicode_New(self->read_len, 127);
    data = (char *)PyUnicode_1BYTE_DATA(ret);
    memcpy(data, self->seq, self->read_len);

    complement_seq(data);

    return ret;
}

PyObject* pyfastx_read_antisense(pyfastx_Read *self, void* closure) {
    char *data;
    PyObject* ret;

    pyfastx_read_get_seq(self);

    ret = PyUnicode_New(self->read_len, 127);
    data = (char *)PyUnicode_1BYTE_DATA(ret);
    memcpy(data, self->seq, self->read_len);

    reverse_complement_seq(data);

    return ret;
}

PyObject* pyfastx_read_description(pyfastx_Read *self, void* closure) {
    Py_ssize_t new_offset;

    if (self->middle->iterating) {
        pyfastx_read_continue_reader(self);
    } else if (!self->desc) {
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
    if (self->middle->iterating) {
        pyfastx_read_continue_reader(self);
    } else if (!self->qual) {
        self->qual = (char *)malloc(self->read_len + 1);
        pyfastx_read_random_reader(self, self->qual, self->qual_offset, self->read_len);
        self->qual[self->read_len] = '\0';
    }

    return Py_BuildValue("s", self->qual);
}

PyObject* pyfastx_read_quali(pyfastx_Read *self, void* closure) {
    int i;
    int phred;

    PyObject *quals;
    PyObject *q;

    if (self->middle->iterating) {
        pyfastx_read_continue_reader(self);
    } else if (!self->qual) {
        self->qual = (char *)malloc(self->read_len + 1);
        pyfastx_read_random_reader(self, self->qual, self->qual_offset, self->read_len);
        self->qual[self->read_len] = '\0';
    }

    phred = self->middle->phred ? self->middle->phred : 33;

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
    {"reverse", (getter)pyfastx_read_reverse, NULL, NULL, NULL},
    {"complement", (getter)pyfastx_read_complement, NULL, NULL, NULL},
    {"antisense", (getter)pyfastx_read_antisense, NULL, NULL, NULL},
    {"qual", (getter)pyfastx_read_qual, NULL, NULL, NULL},
    {"quali", (getter)pyfastx_read_quali, NULL, NULL, NULL},
    {"description", (getter)pyfastx_read_description, NULL, NULL, NULL},
    {NULL}
};

static PyMappingMethods pyfastx_read_as_mapping = {
	.mp_length = (lenfunc)pyfastx_read_length,
};

static PyMemberDef pyfastx_read_members[] = {
    {"id", T_PYSSIZET, offsetof(pyfastx_Read, id), READONLY},
	{"name", T_STRING, offsetof(pyfastx_Read, name), READONLY},
	{NULL}
};

PyTypeObject pyfastx_ReadType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Read",
    .tp_basicsize = sizeof(pyfastx_Read),
    .tp_dealloc = (destructor)pyfastx_read_dealloc,
    .tp_repr = (reprfunc)pyfastx_read_repr,
    .tp_as_mapping = &pyfastx_read_as_mapping,
    .tp_str = (reprfunc)pyfastx_read_str,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_members = pyfastx_read_members,
    .tp_getset = pyfastx_read_getsets,
    .tp_alloc = PyType_GenericAlloc,
    .tp_new = PyType_GenericNew,
};