#include "sequence.h"
#include "structmember.h"


PyObject *pyfastx_sequence_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	pyfastx_Sequence *obj = (pyfastx_Sequence *)type->tp_alloc(type, 0);
	return (PyObject *)obj;
}

int pyfastx_sequence_length(pyfastx_Sequence* self){
	return self->seq_len;
}

PyObject *pyfastx_seqeunce_subscript(pyfastx_Sequence* self, PyObject* item){
	char *seq;
	if (PyIndex_Check(item)) {
		Py_ssize_t i;
		i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		
		if (i == -1 && PyErr_Occurred()){
			return NULL;
		}

		if (i < 0) {
			i += self->seq_len;
		}

		seq = pyfastx_index_get_seq(self->index, self->name, self->offset, self->byte_len, self->start, self->end);
		return Py_BuildValue("C", *(seq+i));
	
	} else if (PySlice_Check(item)) {
		Py_ssize_t slice_start, slice_stop, slice_step, slice_len;

		if(PySlice_Unpack(item, &slice_start, &slice_stop, &slice_step) < 0) {
			return NULL;
		}

		slice_len = PySlice_AdjustIndices(self->seq_len, &slice_start, &slice_stop, slice_step);

		if (slice_len <= 0) {
			Py_RETURN_NONE;
		}

		if(slice_step != 1) {
			Py_RETURN_NONE;
		}

		//create a new sequence
		pyfastx_Sequence *seq = PyObject_New(pyfastx_Sequence, &pyfastx_SequenceType);
		if(!seq){
			return NULL;
		}
		
		seq->name = self->name;
		seq->start = slice_start + 1;
		seq->end = slice_stop;	
		seq->seq_len = slice_stop - slice_start;
		seq->line_len = self->line_len;
		seq->end_len = self->end_len;
		seq->normal = self->normal;
		seq->offset = self->offset;
		seq->byte_len = self->byte_len;
		seq->index = self->index;

		int line_num, tail_num;

		if (self->normal) {
			line_num = seq->seq_len / (seq->line_len - seq->end_len);
			tail_num = seq->seq_len % (seq->line_len - seq->end_len);
			seq->offset = seq->byte_len + seq->start + (seq->start / (seq->line_len - seq->end_len)) * seq->end_len - 1;
			seq->byte_len = line_num * seq->line_len + tail_num;
		}

		printf("%d\n", seq->offset);
		printf("%d\n", seq->byte_len);

		Py_INCREF(seq);
		return (PyObject *)seq;
	} else {
		return NULL;
	}

}

static PyMappingMethods pyfastx_sequence_as_mapping = {
	(lenfunc)pyfastx_sequence_length,
	(binaryfunc)pyfastx_seqeunce_subscript,
	0
};

PyObject *test(pyfastx_Sequence *self, PyObject *args, PyObject *kwargs){
	return Py_BuildValue("i", 100);
}

static PyMethodDef pyfastx_sequence_methods[] = {
	{"test", (PyCFunction)test, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

PyObject *pyfastx_sequence_get_name(pyfastx_Sequence *self, void* closure){
	return Py_BuildValue("s", self->name);
}

PyObject *pyfastx_sequence_get_seq(pyfastx_Sequence *self, void* closure){
	char *seq;
	seq = pyfastx_index_get_seq(self->index, self->name, self->offset, self->byte_len, self->start, self->end);
	return Py_BuildValue("s", seq);
}

static PyGetSetDef pyfastx_sequence_getsets[] = {
	{"name", (getter)pyfastx_sequence_get_name, NULL, NULL, NULL},
	{"seq", (getter)pyfastx_sequence_get_seq, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef pyfastx_sequence_members[] = {
	//{"name", T_STRING, offsetof(pyfastx_Sequence, name), READONLY},
	{"start", T_INT, offsetof(pyfastx_Sequence, start), READONLY},
	{"end", T_INT, offsetof(pyfastx_Sequence, end), READONLY},
	//{"length", T_INT, offsetof(pyfastx_Sequence, seq_len), READONLY},
	{"gc_content", T_FLOAT, offsetof(pyfastx_Sequence, gc_content), READONLY},
	{"composition", T_OBJECT, offsetof(pyfastx_Sequence, composition), READONLY},
	{NULL}
};

PyTypeObject pyfastx_SequenceType = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "Sequence",                     /* tp_name */
    sizeof(pyfastx_Sequence),       /* tp_basicsize */
    0,                              /* tp_itemsize */
    0,                              /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                              /* tp_as_sequence */
    &pyfastx_sequence_as_mapping,   /* tp_as_mapping */
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
    pyfastx_sequence_methods,       /* tp_methods */
    pyfastx_sequence_members,       /* tp_members */
    pyfastx_sequence_getsets,       /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    pyfastx_sequence_new,           /* tp_new */
};
