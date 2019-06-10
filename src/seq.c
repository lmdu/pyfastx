#include "seq.h"

PyObject* pyfastx_sequence_init(SequenceObject* self, PyObject* args, PyObject* kwargs){
    return Py_BuildValue("i", 1);
}

static PyMethodDef sequence_methods[] = {
	{"build_index", (PyCFunction)build_index, METH_VARARGS},
	{"test", (PyCFunction)test, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

PyTypeObject pyfastx_SeqType = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "Sequence",                     /* tp_name */
    sizeof(SequenceObject),         /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)fastx_tp_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                              /* tp_as_sequence */
    0,                              /* tp_as_mapping */
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
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    sequence_methods,               /* tp_methods */
    0,                              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    0,                              /* tp_new */
};