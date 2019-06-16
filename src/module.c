#include "Python.h"
#include "fasta.h"
#include "util.h"
#include "sequence.h"

static PyMethodDef module_methods[] = {
	//{"test", test, METH_VARARGS|METH_KEYWORDS},
	{"clean_seq", clean_seq, METH_VARARGS},
	{"sub_seq", sub_seq, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef module_definition = {
	PyModuleDef_HEAD_INIT,
	"pyfastx",
	"",
	-1,
	module_methods,
};

PyMODINIT_FUNC PyInit_pyfastx(void){
	PyObject *module = PyModule_Create(&module_definition);
	if(!module){
		return NULL;
	}

	if(PyType_Ready(&pyfastx_FastaType) < 0){
		return NULL;
	}
	
	Py_INCREF((PyObject *)&pyfastx_FastaType);
	PyModule_AddObject(module, "Fasta", (PyObject *)&pyfastx_FastaType);

	if(PyType_Ready(&pyfastx_SequenceType) < 0){
		return NULL;
	}

	Py_INCREF((PyObject *)&pyfastx_SequenceType);
	PyModule_AddObject(module, "Sequence", (PyObject *)&pyfastx_SequenceType);

	return module;
}