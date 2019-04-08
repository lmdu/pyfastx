#include "fastx.h"
#include "util.h"

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

    if(PyType_Ready(&pyfastx_FastxType) < 0){
    	return NULL;
    }
    Py_INCREF((PyObject *)&pyfastx_FastxType);
    PyModule_AddObject(module, "fastx", (PyObject *)&pyfastx_FastxType);
    return module;
}