#include "read.h"
#include "structmember.h"

PyObject *pyfastx_read_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
    pyfastx_Read *obj = (pyfastx_Read *)type->tp_alloc(type, 0);
    return (PyObject *)obj;
}

void pyfastx_read_dealloc(pyfastx_Read *self) {
    Py_TYPE(self)->tp_free(self);
}


int pyfastx_read_length(pyfastx_Read *self) {
    return self->read_len;
}

PyObject* pyfastx_read_seq(pyfastx_Read *self, void* closure) {
    if (self->seq == NULL) {
        gzseek(self->ks->f, self->seq_offset, SEEK_SET);
        kstring_t seq = {0, 0, 0};
        if(ks_getuntil(self->ks, '\n', &seq, 0) >= 0) {
           self->seq = seq.s;
        }
    }

    if (self->seq) {
        return Py_BuildValue("s", self->seq);
    }

    return NULL;
}

PyObject* pyfastx_read_qual(pyfastx_Read *self, void* closure) {
    if (self->qual == NULL) {
        gzseek(self->ks->f, self->qual_offset, SEEK_SET);
        kstring_t qual = {0, 0, 0};
        if(ks_getuntil(self->ks, '\n', &qual, 0) >= 0) {
            self->qual = qual.s;
        }
    }

    if (self->qual) {
        return Py_BuildValue("s", self->qual);
    }

    return NULL;
}

static PyGetSetDef pyfastx_read_getsets[] = {
    {"seq", (getter)pyfastx_read_seq, NULL, NULL, NULL},
    {"qual", (getter)pyfastx_read_qual, NULL, NULL, NULL},
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
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                   /* tp_as_sequence */
    &pyfastx_read_as_mapping,                   /* tp_as_mapping */
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
    pyfastx_read_new,              /* tp_new */
};