#ifndef PYFASTX_PYCOMPAT_H
#define PYFASTX_PYCOMPAT_H
#include "Python.h"



//PY_MAJOR_VERSION < 3
#if PY_MAJOR_VERSION < 3
	#define PyExc_ConnectionError PyExc_RuntimeError
	#define PyExc_FileExistsError PyExc_RuntimeError

	#define PyUnicode_AsUTF8 PyString_AsString

	#define PyVarObject_HEAD_INIT(type, size) \
		PyObject_HEAD_INIT(type) size,
	#define PyUnicode_AsUTF8 PyString_AsString
	#define pyfastxSlice_GetIndicesEx(slice, length, start, stop, step, slicelength) \
		PySlice_GetIndicesEx((PySliceObject *)slice, length, start, stop, step, slicelength)
#else
	#define pyfastxSlice_GetIndicesEx PySlice_GetIndicesEx
#endif

#endif