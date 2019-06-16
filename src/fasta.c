#include "fasta.h"
#include "util.h"
#include "structmember.h"

/*calculate fasta attributes including sequence count, length,
composition (ATGCN count) and GC content
*/
void pyfastx_calc_fasta_attrs(pyfastx_Fasta *self){
	//ACGTN nucleotide counts
	int a, c, g, t, n;

	sqlite3_stmt *stmt;
	
	//sequence count
	sqlite3_prepare_v2(self->index->index_db, "SELECT COUNT(*) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	self->seq_counts = sqlite3_column_int(stmt, 0);
	sqlite3_reset(stmt);

	//sequence length
	sqlite3_prepare_v2(self->index->index_db, "SELECT SUM(slen) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	self->seq_length = sqlite3_column_int64(stmt, 0);
	sqlite3_reset(stmt);

	//calculate base counts
	sqlite3_prepare_v2(self->index->index_db, "SELECT SUM(a),SUM(c),SUM(g),SUM(t),SUM(n) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	a = sqlite3_column_int(stmt, 0);
	c = sqlite3_column_int(stmt, 1);
	g = sqlite3_column_int(stmt, 2);
	t = sqlite3_column_int(stmt, 3);
	n = sqlite3_column_int(stmt, 4);
	self->composition = Py_BuildValue("{s:i,s:i,s:i,s:i,s:i}", "A", a, "C", c, "G", g, "T", t, "N", n);
	sqlite3_finalize(stmt);

	//calc GC content
	self->gc_content = (float)(g+c)/(a+c+g+t)*100;
}


PyObject *pyfastx_fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	//fasta file path
	char *file_name;

	//bool value for uppercase sequence
	int uppercase = 1;

	//build index or not
	int build_index = 1;

	//paramters for fasta object construction
	static char* keywords[] = {"file_name", "uppercase", "build_index", NULL};
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s|p|p", keywords, &file_name, &uppercase, &build_index)){
		return NULL;
	}

	//check input sequence file is whether exists
	if(!file_exists(file_name)){
		return PyErr_Format(PyExc_FileExistsError, "input sequence file %s does not exists", file_name);
	}

	//create Fasta class
	pyfastx_Fasta *obj = (pyfastx_Fasta *)type->tp_alloc(type, 0);
	if (!obj){
		return NULL;
	}
	
	//initial sequence file name
	obj->file_name = (char *)malloc(strlen(file_name)+1);
	strcpy(obj->file_name, file_name);

	obj->uppercase = uppercase;

	//create index
	obj->index = pyfastx_init_index(obj->file_name, uppercase);

	//if build_index is True
	if(build_index){
		pyfastx_build_index(obj->index);
		pyfastx_calc_fasta_attrs(obj);
	}
	
	return (PyObject *)obj;
}

void pyfastx_fasta_dealloc(pyfastx_Fasta *self){
	pyfastx_index_free(self->index);
	Py_TYPE(self)->tp_free(self);
}

PyObject *pyfastx_fasta_iter(pyfastx_Fasta *self){
	pyfastx_rewind_index(self->index);
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fasta_next(pyfastx_Fasta *self){
	return pyfastx_get_next_seq(self->index);
}

PyObject *pyfastx_fasta_build_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	pyfastx_build_index(self->index);
	pyfastx_calc_fasta_attrs(self);
	return Py_BuildValue("");
}

PyObject *pyfastx_fasta_rebuild_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	if(file_exists(self->index->index_file)){
		remove(self->index->index_file);
	}
	pyfastx_build_index(self->index);
	pyfastx_calc_fasta_attrs(self);
	return Py_BuildValue("");
}

/*
int fasta_get_item(pyfastx_Fasta *self, PyObject *key){
	if(PyUnicode_Check(key)){
		return 1;
	}
	return 0;
}*/

PyObject *pyfastx_fasta_subscript(pyfastx_Fasta *self, PyObject *item){
	if (PyIndex_Check(item)) {
		Py_ssize_t i;
		i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		
		if (i == -1 && PyErr_Occurred()){
			return NULL;
		}

		if (i < 0) {
			i += self->seq_counts;
		}

		
	}
}

int pyfastx_fasta_length(pyfastx_Fasta *self){
	return self->seq_counts;
}

static PyMemberDef pyfastx_fasta_members[] = {
	{"file_name", T_STRING, offsetof(pyfastx_Fasta, file_name), READONLY},
	{"size", T_LONG, offsetof(pyfastx_Fasta, seq_length), READONLY},
	{"count", T_INT, offsetof(pyfastx_Fasta, seq_counts), READONLY},
	{"gc_content", T_FLOAT, offsetof(pyfastx_Fasta, gc_content), READONLY},
	{"composition", T_OBJECT, offsetof(pyfastx_Fasta, composition), READONLY},
	{NULL}
};

static PyMethodDef pyfastx_fasta_methods[] = {
	{"build_index", (PyCFunction)pyfastx_fasta_build_index, METH_VARARGS},
	{"rebuild_index", (PyCFunction)pyfastx_fasta_rebuild_index, METH_VARARGS},
	//{"get_sub_seq", (PyCFunction)get_sub_seq, METH_VARARGS},
	//{"test", (PyCFunction)test, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

//as a list
static PySequenceMethods seq_methods = {
	0, /*sq_length*/
	(binaryfunc) 0, /*sq_concat*/
	0, /*sq_repeat*/
	0, /*sq_item*/
	0, /*sq_slice */
	0, /*sq_ass_item*/
	0, /*sq_ass_splice*/
	(objobjproc) fasta_get_item, /*sq_contains*/
	(binaryfunc) 0, /*sq_inplace_concat*/
	0,	/*sq_inplace_repeat*/
};

static PyMappingMethods pyfastx_fasta_as_mapping = {
	(lenfunc)pyfastx_fasta_length,
	(binaryfunc)pyfastx_fasta_subscript,
	0,
};

PyTypeObject pyfastx_FastaType = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "Fasta",                        /* tp_name */
    sizeof(pyfastx_Fasta),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_fasta_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    &seq_methods,                   /* tp_as_sequence */
    &pyfastx_fasta_as_mapping,                   /* tp_as_mapping */
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
    (getiterfunc)pyfastx_fasta_iter,     /* tp_iter */
    (iternextfunc)pyfastx_fasta_next,    /* tp_iternext */
    pyfastx_fasta_methods,          /* tp_methods */
    pyfastx_fasta_members,          /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    pyfastx_fasta_new,              /* tp_new */
};
