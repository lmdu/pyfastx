#include <Python.h>
#include "fasta.h"
#include "fastq.h"
#include "fastx.h"
#include "util.h"
#include "read.h"
#include "sequence.h"
#include "fakeys.h"
#include "fqkeys.h"
#include "version.h"
#include "sqlite3.h"
#include "zlib.h"

PyObject *pyfastx_version(PyObject *self, PyObject *args, PyObject *kwargs)	{
	int debug = 0;

	static char* keywords[] = {"debug", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|i", keywords, &debug)) {
		return NULL;
	}

	if (debug) {
		return PyUnicode_FromFormat("pyfastx: %s; zlib: %s; sqlite: %s; zran: %s", PYFASTX_VERSION, ZLIB_VERSION, SQLITE_VERSION, ZRAN_VERSION);
	}

	return Py_BuildValue("s", PYFASTX_VERSION);
}

PyObject *pyfastx_gzip_check(PyObject *self, PyObject *args) {
	PyObject *file_name;

	if (!PyArg_ParseTuple(args, "O", &file_name)) {
		return NULL;
	}

	if (is_gzip_format(file_name)) {
		Py_RETURN_TRUE;
	}

	Py_RETURN_FALSE;
}

PyObject *pyfastx_reverse_complement(PyObject *self, PyObject *args) {
	const char *s;

	PyObject *seq_obj;
	PyObject *rc_obj;

	if (!PyArg_ParseTuple(args, "O", &seq_obj)) {
		return NULL;
	}

	s = PyUnicode_AsUTF8(seq_obj);
	rc_obj = PyUnicode_FromString(s);
	s = PyUnicode_AsUTF8(rc_obj);
	reverse_complement_seq(s);
	return rc_obj;
}

static PyMethodDef module_methods[] = {
	{"version", (PyCFunction)pyfastx_version, METH_VARARGS | METH_KEYWORDS, NULL},
	{"gzip_check", (PyCFunction)pyfastx_gzip_check, METH_VARARGS, NULL},
	{"reverse_complement", (PyCFunction)pyfastx_reverse_complement, METH_VARARGS, NULL},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef module_pyfastx = {
	PyModuleDef_HEAD_INIT,
	"pyfastx",
	"A python C extension for parsing fasta and fastq file",
	-1,
	module_methods,
};


static PyObject* pyfastx_module_init(void){
	PyObject *module = PyModule_Create(&module_pyfastx);

	if (module == NULL){
		return NULL;
	}

	if(PyType_Ready(&pyfastx_FastaType) < 0){
		return NULL;
	}
	Py_INCREF(&pyfastx_FastaType);
	PyModule_AddObject(module, "Fasta", (PyObject *)&pyfastx_FastaType);

	if(PyType_Ready(&pyfastx_FastqType) < 0){
		return NULL;
	}
	Py_INCREF(&pyfastx_FastqType);
	PyModule_AddObject(module, "Fastq", (PyObject *)&pyfastx_FastqType);

	if(PyType_Ready(&pyfastx_FastxType) < 0){
		return NULL;
	}
	Py_INCREF(&pyfastx_FastxType);
	PyModule_AddObject(module, "Fastx", (PyObject *)&pyfastx_FastxType);

	if(PyType_Ready(&pyfastx_SequenceType) < 0){
		return NULL;
	}
	Py_INCREF(&pyfastx_SequenceType);
	PyModule_AddObject(module, "Sequence", (PyObject *)&pyfastx_SequenceType);
	
	if(PyType_Ready(&pyfastx_ReadType) < 0){
		return NULL;
	}
	Py_INCREF(&pyfastx_ReadType);
	PyModule_AddObject(module, "Read", (PyObject *)&pyfastx_ReadType);

	if(PyType_Ready(&pyfastx_FastaKeysType) < 0){
		return NULL;
	}
	Py_INCREF(&pyfastx_FastaKeysType);
	PyModule_AddObject(module, "FastaKeys", (PyObject *)&pyfastx_FastaKeysType);

	if (PyType_Ready(&pyfastx_FastqKeysType) < 0) {
		return NULL;
	}
	Py_INCREF(&pyfastx_FastqKeysType);
	PyModule_AddObject(module, "FastqKeys", (PyObject *)&pyfastx_FastqKeysType);

	PyModule_AddStringConstant(module, "__version__", PYFASTX_VERSION);

	if (!PyErr_Occurred()) {
		return module;
	} else {
		Py_XDECREF(module);
		return NULL;
	}
}

PyMODINIT_FUNC PyInit_pyfastx() {
	return pyfastx_module_init();
}
