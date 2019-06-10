#ifndef PYFASTX_SEQ_H
#define PYFASTX_SEQ_H
#include <Python.h>

//make sequence class
typedef struct {
	PyObject_HEAD
	
	//sequence name
	char* name;

	char* seq;

	//a dict
	PyObject* composition;

} SequenceObject;

extern PyTypeObject pyfastx_SeqType;

#endif