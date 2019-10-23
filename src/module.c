#include "Python.h"
#include "fasta.h"
#include "fastq.h"
#include "util.h"
#include "read.h"
#include "sequence.h"
#include "identifier.h"
#include "version.h"

PyObject *pyfastx_version(PyObject *self, PyObject *args){
	return Py_BuildValue("s", PYFASTX_VERSION);
}

PyObject *pyfastx_gzip_check(PyObject *self, PyObject *args) {
	char *file_name;

	if (!PyArg_ParseTuple(args, "s", &file_name)) {
		return NULL;
	}

	if (is_gzip_format(file_name)) {
		Py_RETURN_TRUE;
	}

	Py_RETURN_FALSE;
}

static PyMethodDef module_methods[] = {
	//{"test", (PyCFunction)pyfastx_test, METH_VARARGS | METH_KEYWORDS, NULL},
	//{"clean_seq", clean_seq, METH_VARARGS, NULL},
	//{"sub_seq", sub_seq, METH_VARARGS, NULL},
	{"version", pyfastx_version, METH_VARARGS, NULL},
	{"gzip_check", pyfastx_gzip_check, METH_VARARGS, NULL},
	{NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef module_pyfastx = {
		PyModuleDef_HEAD_INIT,
		"pyfastx",
		"A python C extension for parsing fasta and fastq file",
		-1,
		module_methods,
	};
#endif

static PyObject* pyfastx_module_init(void){
	PyObject *module;

#if PY_MAJOR_VERSION >= 3
	module = PyModule_Create(&module_pyfastx);

#else
	module = Py_InitModule3("pyfastx", module_methods,
		"A python C extension for parsing fasta and fastq file"
	);
#endif

	if (module == NULL){
		return NULL;
	}

	if(PyType_Ready(&pyfastx_FastaType) < 0){
		return NULL;
	}
	Py_INCREF((PyObject *)&pyfastx_FastaType);
	PyModule_AddObject(module, "Fasta", (PyObject *)&pyfastx_FastaType);

	if(PyType_Ready(&pyfastx_FastqType) < 0){
		return NULL;
	}
	Py_INCREF((PyObject *)&pyfastx_FastqType);
	PyModule_AddObject(module, "Fastq", (PyObject *)&pyfastx_FastqType);

	if(PyType_Ready(&pyfastx_SequenceType) < 0){
		return NULL;
	}
	Py_INCREF((PyObject *)&pyfastx_SequenceType);
	PyModule_AddObject(module, "Sequence", (PyObject *)&pyfastx_SequenceType);

	
	if(PyType_Ready(&pyfastx_ReadType) < 0){
		return NULL;
	}
	Py_INCREF((PyObject *)&pyfastx_ReadType);
	PyModule_AddObject(module, "Read", (PyObject *)&pyfastx_ReadType);

	if(PyType_Ready(&pyfastx_IdentifierType) < 0){
		return NULL;
	}
	Py_INCREF((PyObject *)&pyfastx_IdentifierType);
	PyModule_AddObject(module, "Identifier", (PyObject *)&pyfastx_IdentifierType);

	return module;
}

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_pyfastx(void) {
	return pyfastx_module_init();
}
#else
PyMODINIT_FUNC initpyfastx(void) {
	pyfastx_module_init();
}
#endif