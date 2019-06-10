#include "fasta.h"
#include "util.h"
#include "structmember.h"
KSEQ_INIT(gzFile, gzread, gzrewind)


PyObject *fasta_tp_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	//fasta file path
	char *file_name;

	//bool value for upper sequence
	int uppercase = 0;
	
	//paramters for fasta object construction
	static char* keywords[] = {"file_name", "uppercase", NULL};
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s|p", keywords, &file_name, &uppercase)){
		return NULL;
	}

	//check input sequence file is whether exists
	if(!file_exists(file_name)){
		return PyErr_Format(PyExc_FileExistsError, "input sequence file %s does not exists", file_name);
	}

	//create Fastx class
	FastaObject *obj = (FastaObject *)type->tp_alloc(type, 0);
	if (!obj){
		return NULL;
	}
	
	//initial sequence file name
	obj->file_name = (char *)malloc(strlen(file_name)+1);
	strcpy(obj->file_name, file_name);

	obj->uppercase = uppercase;

	//initial kseqs
	obj->gzfp = gzopen(obj->file_name, "rb");
	obj->kseqs = kseq_init(obj->gzfp);

	//open file and check file is gzip
	obj->fd = fopen(obj->file_name, "rb");
	obj->is_gzip = is_gzip(obj->fd);

	//create index file path
	obj->index_file = (char *)malloc(strlen(file_name) + 4);
	strcpy(obj->index_file, file_name);
	strcat(obj->index_file, ".db");

	//if index file exists, connect to it
	obj->db = NULL;
	obj->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
	if (file_exists(obj->index_file)) {
		if(sqlite3_open(obj->index_file, &obj->db) != SQLITE_OK){
			PyErr_SetString(PyExc_ConnectionError, sqlite3_errmsg(obj->db));
			return NULL;
		}
		load_gzip_index(obj);
	}
	return (PyObject *)obj;
}

void fastx_tp_dealloc(FastxObject *self){
	kseq_destroy(self->kseqs);
	//zran_free(self->gzip_index);
	gzclose(self->gzfp);
	fclose(self->fd);
	Py_TYPE(self)->tp_free(self);
}

PyObject *fastx_tp_iter(FastxObject *self){
	kseq_rewind(self->kseqs);
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *fastx_tp_next(FastxObject *self){
	if(kseq_read(self->kseqs) >= 0){
		if(self->uppercase){
			upper_string(self->kseqs->seq.s);
		}
		return Py_BuildValue("(ss)", self->kseqs->name.s, self->kseqs->seq.s);
	}
	return NULL;
}

int fastx_get_item(FastxObject *self, PyObject *key){
	if(PyUnicode_Check(key)){
		return 1;
	}
	return 0;
}

PyObject *fastx_get_key(FastxObject *self, PyObject *key){
	if(PyLong_Check(key)){
		int index = PyLong_AsLong(key);
		return Py_BuildValue("i", index);
	}

	else if(PyUnicode_Check(key)){
		char *name = PyUnicode_AsUTF8(key);
		return Py_BuildValue("s", name);
	}

	else {
		PyErr_SetObject(PyExc_KeyError, key);
		return NULL;
	}
}

int fastx_get_len(FastxObject *self){
	return 123; 
}

int fastx_get_val(FastxObject *self, PyObject *key, PyObject *val){ 
	char *name = PyUnicode_AsUTF8(key);
	return 1;
}



/*
static PyMemberDef fastx_members[] = {
	{"file_path", T_STRING, offsetof(FastxObject, file_path), 0, "file path"},
	{NULL}
};*/

static PyMethodDef fastx_methods[] = {
	{"build_index", (PyCFunction)build_index, METH_VARARGS},
	{"get_sub_seq", (PyCFunction)get_sub_seq, METH_VARARGS|METH_KEYWORDS},
	{"test", (PyCFunction)test, METH_VARARGS},
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
	(objobjproc) fastx_get_item, /*sq_contains*/
	(binaryfunc) 0, /*sq_inplace_concat*/
	0,	/*sq_inplace_repeat*/
};

static PyMappingMethods map_methods = {
	(lenfunc)fastx_get_len,
	(binaryfunc)fastx_get_key,
	(objobjargproc)fastx_get_val,
};

PyTypeObject pyfastx_FastaType = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "Fasta",                        /* tp_name */
    sizeof(FastxObject),             /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)fastx_tp_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    &seq_methods,                   /* tp_as_sequence */
    &map_methods,                   /* tp_as_mapping */
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
    (getiterfunc)fastx_tp_iter,     /* tp_iter */
    (iternextfunc)fastx_tp_next,    /* tp_iternext */
    fastx_methods,                  /* tp_methods */
    0,                              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    fastx_tp_new,                   /* tp_new */
};
