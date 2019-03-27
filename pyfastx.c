#include <zlib.h>
#include <Python.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

kseq_t *SEQ;


static char upper(char c){
	if(c>='a' && c<='z'){
		return 'A'+c-'a';
	}else{
		return c;
	}
}

static int upper_string(char *str){
	int i;
	for(i=0; str[i]; i++){
		str[i] = upper(str[i]);
	}
	return i;
}

static PyObject *open_fasta(PyObject *self, PyObject *args){
	char *fasta_path;
	if (!PyArg_ParseTuple(args, "s", &fasta_path)){
		return NULL;
	}
	gzFile fp;
	fp = gzopen(fasta_path, "rb");
	SEQ = kseq_init(fp);
	return Py_BuildValue("i", 1);
}

static PyObject *close_fasta(PyObject *self, PyObject *args){
	if(SEQ){
		kseq_destroy(SEQ);
	}
	return Py_BuildValue("i", 1);
}

static PyObject *clean_seq(PyObject *self, PyObject *args){
	char *seq;
	if (!PyArg_ParseTuple(args, "s", &seq)){
		return NULL;
	}
	int i;
	int j = 0;
	for(i=0; seq[i]; i++){
		if(!isspace(seq[i])){
			seq[j++] = upper(seq[i]);
		}
	}
	seq[j] = '\0';
	return Py_BuildValue("s", seq);
}

static PyObject *sub_seq(PyObject *self, PyObject *args){
	char *seq;
	int start;
	int end;

	if (!PyArg_ParseTuple(args, "sii", &seq, &start, &end)){
		return NULL;
	}
	int i;
	int j = 0;
	int real_pos = 0;
	int flag;
	for(i=0; seq[i]; i++){
		flag = isspace(seq[i]);

		if(!flag){
			real_pos++;
		}

		if(real_pos > end){
			break;
		}

		if(real_pos >= start){
			if(!flag){
				seq[j++] = upper(seq[i]);
			}
		}
	}
	seq[j] = '\0';
	return Py_BuildValue("s", seq);
}

static PyObject *iter_seq(PyObject *self, PyObject *args){
	int l;
	while((l=kseq_read(SEQ))>=0){
		upper_string(SEQ->seq.s);
		return Py_BuildValue("(ss)", SEQ->name.s, SEQ->seq.s);
	}
	return Py_BuildValue("");
}

static PyObject *build_index(PyObject *self, PyObject *args){
	char *fasta_path;
	if (!PyArg_ParseTuple(args, "s", &fasta_path)){
		return NULL;
	}
	int position = 0;
	int start = 0;
	int seqlen = 0;
	int gc = 0;
	int ns = 0;
	int c;
	kseq_t *seq;
	kstream_t *ks;
	gzFile fp;
	PyObject *result = PyList_New(0);
	PyObject *tmp;
	
	fp = gzopen(fasta_path, "rb");
	seq = kseq_init(fp);
	ks = seq->f;
	while((c=ks_getc(ks))!=-1){
		position++;
		if(c == 62){
			if(start){
				tmp = Py_BuildValue("(siiiii)", seq->name.s, start, position-start-1, seqlen, gc, ns);
				PyList_Append(result, tmp);
				Py_DECREF(tmp);
			}
			position += ks_getuntil(ks, 0, &seq->name, &c);
			position++;
			while(c != 10){
				c = ks_getc(ks);
				position++;
			}
			start = position;
			seqlen = 0;
			gc = 0;
			ns = 0;
		}else{
			if(c != 10 && c != 13){
				seqlen++;
				c = toupper(c);
				if(c == 71 || c == 67){
					gc++;
				}else if(c == 78){
					ns++;
				}
			}
		}
	}
	tmp = Py_BuildValue("(siiiii)", seq->name.s, start, position-start, seqlen, gc, ns);
	PyList_Append(result, tmp);
	Py_DECREF(tmp);
	return result;
}

//make sequence iterator
typedef struct {
	PyObject_HEAD
	gzFile fp;
	kseq_t *sequence;
} FastaState;

static PyObject *fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	char *fasta_path;
	if (!PyArg_ParseTuple(args, "s", &fasta_path)){
		return NULL;
	}
	gzFile fp;
	kseq_t *sequence;

	fp = gzopen(fasta_path, "rb");
	sequence = kseq_init(fp);
	
	FastaState *fstate = (FastaState *)type->tp_alloc(type, 0);
	if (!fstate){
		return NULL;
	}
	fstate->fp = fp;
	fstate->sequence = sequence;

	return (PyObject *)fstate;
}

static void fasta_dealloc(FastaState *fstate){
	kseq_destroy(fstate->sequence);
	gzclose(fstate->fp);
	Py_TYPE(fstate)->tp_free(fstate);
}

static PyObject *fasta_next(FastaState *fstate){
	int l;
	//if((l=kseq_read(fstate->sequence))>=0){
	//	upper_string(fstate->sequence->seq.s);
	//	return Py_BuildValue("(ss)", fstate->sequence->name.s, fstate->sequence->seq.s);
	//}
	if((l=kseq_read(fstate->sequence))>=0){
		upper_string(fstate->sequence->seq.s);
		return Py_BuildValue("(ss)", fstate->sequence->name.s, fstate->sequence->seq.s);
	}

	return NULL;
}

PyTypeObject PyFasta_Type = {
	PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "fasta",                        /* tp_name */
    sizeof(FastaState),             /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)fasta_dealloc,      /* tp_dealloc */
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
    PyObject_SelfIter,              /* tp_iter */
    (iternextfunc)fasta_next,       /* tp_iternext */
    0,                              /* tp_methods */
    0,                              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    fasta_new,                      /* tp_new */
};

static PyMethodDef kseq_methods[] = {
	{"build_index", build_index, METH_VARARGS},
	{"open_fasta", open_fasta, METH_VARARGS},
	{"close_fasta", close_fasta, METH_VARARGS},
	{"iter_seq", iter_seq, METH_VARARGS},
	{"clean_seq", clean_seq, METH_VARARGS},
	{"sub_seq", sub_seq, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef kseq_definition = {
	PyModuleDef_HEAD_INIT,
	"pyfastx",
	"",
	-1,
	kseq_methods,
};

PyMODINIT_FUNC PyInit_pyfastx(void){
    PyObject *module = PyModule_Create(&kseq_definition);
    if(!module){
    	return NULL;
    }

    if(PyType_Ready(&PyFasta_Type) < 0){
    	return NULL;
    }
    Py_INCREF((PyObject *)&PyFasta_Type);
    PyModule_AddObject(module, "fasta", (PyObject *)&PyFasta_Type);
    return module;
}
