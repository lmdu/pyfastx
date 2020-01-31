#include "read.h"
#include "util.h"
#include "structmember.h"

void pyfastx_read_dealloc(pyfastx_Read *self) {
    if (self->seq) {
        free(self->seq);
    }

    if (self->qual) {
        free(self->qual);
    }

    Py_TYPE(self)->tp_free((PyObject *)self);
}


int pyfastx_read_length(pyfastx_Read *self) {
    return self->read_len;
}

PyObject* pyfastx_read_raw(pyfastx_Read *self, void* closure) {
    int64_t new_offset;
    int64_t new_bytelen;
    char *buff;

    new_offset = self->seq_offset - self->desc_len - 1;
    new_bytelen = self->qual_offset + self->read_len - new_offset + 1;

    buff = (char *)malloc(new_bytelen + 2);

    if (self->gzip_format) {
        zran_seek(self->gzip_index, new_offset, SEEK_SET, NULL);
        zran_read(self->gzip_index, buff, new_bytelen);
    } else {
        gzseek(self->gzfd, new_offset, SEEK_SET);
        gzread(self->gzfd, buff, new_bytelen);
    }

    if (buff[new_bytelen-1] == '\r') {
        buff[new_bytelen] = '\n';
        buff[new_bytelen+1] = '\0';
    } else {
        buff[new_bytelen] = '\0';
    }

    return Py_BuildValue("s", buff);
}

PyObject* pyfastx_read_seq(pyfastx_Read *self, void* closure) {
    if (self->seq == NULL) {
        self->seq = (char *)malloc(self->read_len + 1);

        if (self->gzip_format) {
            zran_seek(self->gzip_index, self->seq_offset, SEEK_SET, NULL);
            zran_read(self->gzip_index, self->seq, self->read_len);
        } else {
            gzseek(self->gzfd, self->seq_offset, SEEK_SET);
            gzread(self->gzfd, self->seq, self->read_len);
        }
        self->seq[self->read_len] = '\0';
    }

    if (self->seq) {
        return Py_BuildValue("s", self->seq);
    }

    Py_RETURN_NONE;
}

PyObject* pyfastx_read_description(pyfastx_Read *self, void* closure) {
    int64_t new_offset;
    char *buff;
    
    new_offset = self->seq_offset - self->desc_len - 1;
    buff = (char *)malloc(self->desc_len + 1);

    if (self->gzip_format) {
        zran_seek(self->gzip_index, new_offset, SEEK_SET, NULL);
        zran_read(self->gzip_index, buff, self->desc_len);
    } else {
        gzseek(self->gzfd, new_offset, SEEK_SET);
        gzread(self->gzfd, buff, self->desc_len);
    }

    if (buff[self->desc_len-1] == '\r') {
        buff[self->desc_len-1] = '\0';
    } else {
        buff[self->desc_len] = '\0';
    }

    return Py_BuildValue("s", buff);
}

PyObject* pyfastx_read_qual(pyfastx_Read *self, void* closure) {
    if (self->qual == NULL) {
        self->qual = (char *)malloc(self->read_len + 1);
        if (self->gzip_format) {
            zran_seek(self->gzip_index, self->qual_offset, SEEK_SET, NULL);
            zran_read(self->gzip_index, self->qual, self->read_len);
        } else {
            gzseek(self->gzfd, self->qual_offset, SEEK_SET);
            gzread(self->gzfd, self->qual, self->read_len);
        }
        self->qual[self->read_len] = '\0';
    }

    if (self->qual) {
        return Py_BuildValue("s", self->qual);
    }

    Py_RETURN_NONE;
}

PyObject* pyfastx_read_quali(pyfastx_Read *self, void* closure) {
    int phred;

    if (self->qual == NULL) {
        pyfastx_read_qual(self, NULL);
    }

    phred = self->phred ? self->phred : 33;

    if (self->qual != NULL) {
        PyObject *quals = PyList_New(0);
        int i;
        for (i = 0; i < self->read_len; i++) {
            PyObject *q = Py_BuildValue("i", self->qual[i] - phred);
            PyList_Append(quals, q);
            Py_DECREF(q);
        }

        return quals;
    }

    return NULL;
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
    {"id", T_LONG, offsetof(pyfastx_Read, id), READONLY},
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