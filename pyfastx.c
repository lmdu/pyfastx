#include <zlib.h>
#include <sqlite3.h>
#include <Python.h>
#include "kseq.h"
#include "zran.h"
//#include "structmember.h"

KSEQ_INIT(gzFile, gzread)

//make sequence iterator
typedef struct {
	PyObject_HEAD
	gzFile gzfp;
	kseq_t* kseqs;
	sqlite3* db;
	char* file_path;
} FastxObject;

int file_exists(char *filename){
	FILE *fp;
	if((fp=fopen(filename, "r"))){
		fclose(fp);
		return 1;
	}
	return 0;
}

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

//check input file is whether gzip file
static int is_gzip(FILE *fd){
	unsigned char magic[4] = {0};
	int ret;

	ret = fread(magic, 1, sizeof(magic), fd);
	fseek(fd, 0, SEEK_SET);

	if(ret != sizeof(magic)){
		return 0;
	}
	if(magic[0] != 0x1f || magic[1] != 0x8b || magic[2] != 0x08){
		return 0;
	}
	return 1;
}

static PyObject *build_index(FastxObject *self, PyObject *args, PyObject *kwargs){
	static char* keywords[] = {"rebuild", NULL};
	int rebuild = 0;
	int ret;
	char *index_file;
	sqlite3_stmt *stmt;
	
	/*1: normal fasta sequence with the same length in line
	  0: not normal fasta sequence with different length in line*/
	int seq_normal = 1;

	//1: \n, 2: \r\n
	int line_end = 1;
	
	//length of previous line
	int line_len = 0;

	//length of current line
	int temp_len = 0;

	//number of lines that line_len not equal to temp_len
	int bad_line = 0;

	int position = 0;
	int start = 0;
	int seq_len = 0;
	int gc = 0;
	int ns = 0;
	int c;
	kseq_t *seq;
	kstream_t *ks;
	gzFile fp;

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "|p", keywords, &rebuild)){
		return NULL;
	}

	//create index file path
	index_file = malloc(strlen(self->file_path) + 4);
	strcpy(index_file, self->file_path);
	strcat(index_file, ".db");

	//check index file exists, if exists do not build index
	if(file_exists(index_file)){
		if(rebuild){
			if(remove(index_file)==0){
				printf("removed %s\n", index_file);
			}
		} else {
			return Py_BuildValue("i", 1);
		}
	}

	ret = sqlite3_open(index_file, &self->db);
	if(ret != SQLITE_OK){
		return Py_BuildValue("i", 0);
	}

	//create index database
	const char *create_sql = "\
		CREATE TABLE seq (\
			seqid TEXTPRIMARY KEY,\
			start INTEGER,\
			end INTEGER,\
			seqlen INTEGER,\
			linelen INTEGER,\
			lineend INTEGER,\
			normal INTEGER,\
			gc INTEGER,\
			ns INTEGER\
		);\
		CREATE TABLE gzindex (\
			ID INTEGER PRIMARY KEY,\
			content BLOB\
		);";

	char *errmsg;
	ret = sqlite3_exec(self->db, create_sql, NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "create table error");
	}

	ret = sqlite3_exec(self->db, "PRAGMA synchronous=OFF", NULL, NULL, &errmsg);
	if(ret != SQLITE_OK){
		printf("%s\n", "pragma");
		printf("%s\n", errmsg);
	}

	ret = sqlite3_exec(self->db, "begin", NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "begin error");
	}

	char *insert_sql = "INSERT INTO seq VALUES (?,?,?,?,?,?,?,?,?)";
	
	sqlite3_prepare_v2(self->db, insert_sql, -1, &stmt, NULL);

	fp = gzopen(self->file_path, "rb");
	seq = kseq_init(fp);
	ks = seq->f;
	while((c=ks_getc(ks))!=-1){
		position++;
		
		// c is >
		if(c == 62){
			if(start){

				//end of sequenc and check whether normal fasta
				if(bad_line > 1){
					seq_normal = 0;
				}
				sqlite3_bind_text(stmt, 1, seq->name.s, seq->name.l, NULL);
				sqlite3_bind_int(stmt, 2, start);
				sqlite3_bind_int(stmt, 3, position-start-1);
				sqlite3_bind_int(stmt, 4, seq_len);
				sqlite3_bind_int(stmt, 5, line_len);
				sqlite3_bind_int(stmt, 6, line_end);
				sqlite3_bind_int(stmt, 7, seq_normal);
				sqlite3_bind_int(stmt, 8, gc);
				sqlite3_bind_int(stmt, 9, ns);
				sqlite3_step(stmt);
				sqlite3_reset(stmt);
				printf("%s\n", "yes");
			}
			position += ks_getuntil(ks, 0, &seq->name, &c);
			position++;
			while(c != 10){
				c = ks_getc(ks);
				position++;
			}
			start = position;
			seq_len = 0;
			gc = 0;
			ns = 0;
			temp_len = 0;
			line_len = 0;
			line_end = 1;
			seq_normal = 1;
		}

 		// c is \r
		else if(c == 13){
			temp_len++;

			if(line_end != 2){
				line_end = 2;
			}
		}
		
		// c is \n
		else if(c == 10){
			temp_len++;
			if(line_len){
				if(line_len != temp_len){
					bad_line++;
				}
			} else {
				line_len = temp_len;
			}

		}

		else {
			seq_len++;

			//temp line length
			temp_len++;

			c = toupper(c);
			
			// c is G or C
			if(c == 71 || c == 67){
				gc++;
			}
			// c is N, unkown base
			else if(c == 78){
				ns++;
			}
		}
	}

	//end of sequenc and check whether normal fasta
	if(bad_line > 1){
		seq_normal = 0;
	}

	sqlite3_bind_text(stmt, 1, seq->name.s, seq->name.l, NULL);
	sqlite3_bind_int(stmt, 2, start);
	sqlite3_bind_int(stmt, 3, position-start-1);
	sqlite3_bind_int(stmt, 4, seq_len);
	sqlite3_bind_int(stmt, 5, line_len);
	sqlite3_bind_int(stmt, 6, line_end);
	sqlite3_bind_int(stmt, 7, seq_normal);
	sqlite3_bind_int(stmt, 8, gc);
	sqlite3_bind_int(stmt, 9, ns);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	sqlite3_exec(self->db, "commit;", NULL, NULL, NULL);
	
	return Py_BuildValue("i", 1);
}

static PyObject *test(PyObject *self, PyObject *args){
	/*static char* keywords[] = {"fasta", NULL};
	char* fasta;
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s", keywords, &fasta)){
		return NULL;
	}*/
	char *fasta_path;
	if (!PyArg_ParseTuple(args, "s", &fasta_path)){
		return NULL;
	}

	FILE *fp = fopen(fasta_path, "rb");

	int g = is_gzip(fp);
	printf("is gzip: %d\n", g);
	return Py_BuildValue("");
	printf("%s\n", fasta_path);
	
	FILE *fd = fopen("idex.idx", "wb");

	zran_index_t *index = (zran_index_t *)malloc(sizeof(zran_index_t));

	zran_init(index, fp, 0, 0, 0, 0);
	zran_build_index(index, 0, 0);
	zran_export_index(index, fd);

	sqlite3 *db = NULL;
	sqlite3_open("index.db", &db);
	sqlite3_exec(db, "CREATE TABLE seqidx (content blob);", NULL, NULL, NULL);
	sqlite3_stmt *stmt;
	sqlite3_prepare_v2(db, "INSERT INTO seqidx VALUES (?)", -1, &stmt, NULL);
	sqlite3_bind_blob(stmt, 1, index, sizeof(index), NULL);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
	sqlite3_close(db);
	
	fclose(fp);
	fclose(fd);
	
	return Py_BuildValue("i", 1);

	/*
	if(c){
		return Py_BuildValue("i", a+b);
	} else {
		return Py_BuildValue("i", a-b);
	}*/
}

static PyObject *fastx_tp_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	char *file_path;
	if (!PyArg_ParseTuple(args, "s", &file_path)){
		return NULL;
	}
	gzFile gzfp;
	kseq_t *kseqs;

	gzfp = gzopen(file_path, "rb");
	kseqs = kseq_init(gzfp);
	
	FastxObject *self = (FastxObject *)type->tp_alloc(type, 0);
	if (!self){
		return NULL;
	}
	self->file_path = (char *)malloc(strlen(file_path)+1);
	strcpy(self->file_path, file_path);
	self->gzfp = gzfp;
	self->kseqs = kseqs;
	self->db = NULL;

	printf("%s\n", self->file_path);

	return (PyObject *)self;
}

static void fastx_tp_dealloc(FastxObject *slef){
	kseq_destroy(slef->kseqs);
	gzclose(slef->gzfp);
	Py_TYPE(slef)->tp_free(slef);
}

static PyObject *fastx_tp_next(FastxObject *slef){
	int l;
	if((l=kseq_read(slef->kseqs))>=0){
		upper_string(slef->kseqs->seq.s);
		return Py_BuildValue("(ss)", slef->kseqs->name.s, slef->kseqs->seq.s);
	}

	return NULL;
}

/*
static PyMemberDef fastx_members[] = {
	{"file_path", T_STRING, offsetof(FastxObject, file_path), 0, "file path"},
	{NULL}
};*/

static PyMethodDef fastx_methods[] = {
	{"build_index", (PyCFunction)build_index, METH_VARARGS},
	{"clean_seq", clean_seq, METH_VARARGS},
	{"sub_seq", sub_seq, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

static PyMethodDef module_methods[] = {
	//{"test", test, METH_VARARGS|METH_KEYWORDS},
	{"test", test, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

PyTypeObject PyFastx_Type = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "fastx",                        /* tp_name */
    sizeof(FastxObject),             /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)fastx_tp_dealloc,   /* tp_dealloc */
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


static struct PyModuleDef kseq_definition = {
	PyModuleDef_HEAD_INIT,
	"pyfastx",
	"",
	-1,
	module_methods,
};

PyMODINIT_FUNC PyInit_pyfastx(void){
    PyObject *module = PyModule_Create(&kseq_definition);
    if(!module){
    	return NULL;
    }

    if(PyType_Ready(&PyFastx_Type) < 0){
    	return NULL;
    }
    Py_INCREF((PyObject *)&PyFastx_Type);
    PyModule_AddObject(module, "fastx", (PyObject *)&PyFastx_Type);
    return module;
}
