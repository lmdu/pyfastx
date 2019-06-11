#include "fasta.h"
#include "util.h"
#include "structmember.h"
KSEQ_INIT(gzFile, gzread, gzrewind)


void _pyfastx_build_gzip_index(pyfastx_Fasta *self){
	sqlite3_stmt *stmt;

	zran_init(self->gzip_index, self->fd, 0, 0, 0, ZRAN_AUTO_BUILD);
	zran_build_index(self->gzip_index, 0, 0);

	//create temp gzip index file
	char *temp_index = (char *)malloc(strlen(self->file_name) + 5);
	strcpy(temp_index, self->file_name);
	strcat(temp_index, ".tmp");
	FILE* fd = fopen(temp_index, "wb");
	zran_export_index(self->gzip_index, fd);
	
	long fsize = ftell(fd);
	fseek(fd, 0, SEEK_SET);
	char *buff = (char *)malloc(fsize + 1);

	if(fread(buff, 1, fsize, fd) != 0){
		PyErr_SetString(PyExc_RuntimeError, "Error occured when reading gzip index from temp file.");
		return;
	}
	fclose(fd);
	remove(temp_index);

	sqlite3_prepare_v2(self->index_db, "INSERT INTO gzindex VALUES (NULL, ?)", -1, &stmt, NULL);
	sqlite3_bind_blob(stmt, 1, buff, strlen(buff), NULL);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
}

void _pyfastx_load_gzip_index(pyfastx_Fasta *self){
	sqlite3_stmt *stmt;
	zran_init(self->gzip_index, self->fd, 0, 0, 0, ZRAN_AUTO_BUILD);
	char *temp_index = (char *)malloc(strlen(self->file_name) + 5);
	strcpy(temp_index, self->file_name);
	strcat(temp_index, ".tmp");
	FILE* fd = fopen(temp_index, "wb+");

	sqlite3_prepare_v2(self->index_db, "SELECT content FROM gzindex;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	const char *buff = sqlite3_column_blob(stmt, 0);
	fwrite(buff, 1, strlen(buff), fd);
	fseek(fd, 0, SEEK_SET);

	zran_import_index(self->gzip_index, fd);

	fclose(fd);
	remove(temp_index);
}

/*calculate fasta attributes including sequence count, length,
composition (ATGCN count) and GC content
*/
void _pyfastx_calc_fasta_attrs(pyfastx_Fasta *self){
	int a_counts;
	int c_counts;
	int g_counts;
	int t_counts;
	int n_counts;

	sqlite3_stmt *stmt;
	
	//sequence count
	sqlite3_prepare_v2(self->index_db, "SELECT COUNT(*) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	self->seq_counts = sqlite3_column_int(stmt, 0);
	sqlite3_reset(stmt);

	//sequence length
	sqlite3_prepare_v2(self->index_db, "SELECT SUM(slen) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	self->seq_length = sqlite3_column_int64(stmt, 0);
	sqlite3_reset(stmt);

	//calculate base counts
	sqlite3_prepare_v2(self->index_db, "SELECT SUM(a),SUM(c),SUM(g),SUM(t),SUM(n) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	a_counts = sqlite3_column_int(stmt, 0);
	c_counts = sqlite3_column_int(stmt, 1);
	g_counts = sqlite3_column_int(stmt, 2);
	t_counts = sqlite3_column_int(stmt, 3);
	n_counts = sqlite3_column_int(stmt, 4);
	self->composition = Py_BuildValue("{s:i,s:i,s:i,s:i,s:i}", "A", a_counts, "C", c_counts, "G", g_counts, "T", t_counts, "N", n_counts);
	sqlite3_finalize(stmt);

	//calc GC content
	self->gc_content = (float)(g_counts+c_counts)/(a_counts+c_counts+g_counts+t_counts)*100;
}

void _pyfastx_build_index(pyfastx_Fasta *self){
	//time
	time_t start_time, current_time;

	start_time = time(NULL);

	// seqlite3 return value
	int ret;
	
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	
	// 1: normal fasta sequence with the same length in line
	// 0: not normal fasta sequence with different length in line
	int seq_normal = 1;

	//1: \n, 2: \r\n
	int line_end = 1;
	
	//length of previous line
	int line_len = 0;

	//length of current line
	int temp_len = 0;

	//number of lines that line_len not equal to temp_len
	int bad_line = 0;

	//current read position
	int position = 0;
	
	// start position
	int start = 0;

	//sequence length
	int seq_len = 0;

	//number of bases
	int g_count = 0;
	int c_count = 0;
	int a_count = 0;
	int t_count = 0;
	int n_count = 0;

	//current read base char
	int c;

	//reading file for kseq
	kstream_t *ks;


	ret = sqlite3_open(self->index_file, &self->index_db);
	if(ret != SQLITE_OK){
		PyErr_SetString(PyExc_ConnectionError, sqlite3_errmsg(self->index_db));
		return;
	}

	//create index database
	const char *create_sql = " \
		CREATE TABLE seq ( \
			ID INTEGER PRIMARY KEY, --seq identifier\n \
			seqid TEXT, --seq name\n \
			offset INTEGER, --seq offset start\n \
			blen INTEGER, --seq byte length\n \
			slen INTEGER, --seq length\n \
			llen INTEGER, --line lenght\n \
			elen INTEGER, --end length\n \
			norm INTEGER, --line with same length or not\n \
			a INTEGER, --A base counts\n \
			c INTEGER, --C base counts\n \
			g INTEGER, --G base counts\n \
			t INTEGER, --T base counts\n \
			n INTEGER --unknown base counts\n \
		);\
		CREATE TABLE gzindex ( \
			ID INTEGER PRIMARY KEY, \
			content BLOB \
		);";

	ret = sqlite3_exec(self->index_db, create_sql, NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "create table error");
	}

	ret = sqlite3_exec(self->index_db, "PRAGMA synchronous=OFF;", NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "pragma");
	}

	ret = sqlite3_exec(self->index_db, "BEGIN;", NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "begin error");
	}

	const char *insert_sql = "INSERT INTO seq VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)";
	
	sqlite3_prepare_v2(self->index_db, insert_sql, -1, &stmt, NULL);

	current_time = time(NULL);
	printf("time: %ld\n", current_time - start_time);
	start_time = current_time;

	ks = self->kseqs->f;
	while((c=ks_getc(ks))!=-1){
		position++;
		
		// c is >
		if(c == 62){
			if(start){

				//end of sequenc and check whether normal fasta
				if(bad_line > 1){
					seq_normal = 0;
				}

				current_time = time(NULL);
				printf("%ld\n", current_time-start_time);
				start_time = current_time;

				sqlite3_bind_null(stmt, 1);
				sqlite3_bind_text(stmt, 2, self->kseqs->name.s, -1, NULL);
				sqlite3_bind_int(stmt, 3, start);
				sqlite3_bind_int(stmt, 4, position-start-1);
				sqlite3_bind_int(stmt, 5, seq_len);
				sqlite3_bind_int(stmt, 6, line_len);
				sqlite3_bind_int(stmt, 7, line_end);
				sqlite3_bind_int(stmt, 8, seq_normal);
				sqlite3_bind_int(stmt, 9, a_count);
				sqlite3_bind_int(stmt, 10, c_count);
				sqlite3_bind_int(stmt, 11, g_count);
				sqlite3_bind_int(stmt, 12, t_count);
				sqlite3_bind_int(stmt, 13, n_count);
				sqlite3_step(stmt);
				sqlite3_reset(stmt);

				current_time = time(NULL);
				printf("%ld\n", current_time-start_time);
				current_time = start_time;
			}
			position += ks_getuntil(ks, 0, &self->kseqs->name, &c);
			position++;
			while(c != 10){
				c = ks_getc(ks);
				position++;
			}
			start = position;
			seq_len = 0;
			g_count = 0;
			c_count = 0;
			a_count = 0;
			t_count = 0;
			n_count = 0;
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
			
			//calculate base counts in sequence
			switch(c){
				case 65:
					a_count++; break;
				
				case 67:
					c_count++; break;
				
				case 71:
					g_count++; break;
				
				case 84:
					t_count++; break;

				default:
					n_count++;
			}
		}
	}

	//reset read position of sequence file
	kseq_rewind(self->kseqs);

	//end of sequenc and check whether normal fasta
	if(bad_line > 1){
		seq_normal = 0;
	}

	sqlite3_bind_null(stmt, 1);
	sqlite3_bind_text(stmt, 2, self->kseqs->name.s, -1, NULL);
	sqlite3_bind_int(stmt, 3, start);
	sqlite3_bind_int(stmt, 4, position-start-1);
	sqlite3_bind_int(stmt, 5, seq_len);
	sqlite3_bind_int(stmt, 6, line_len);
	sqlite3_bind_int(stmt, 7, line_end);
	sqlite3_bind_int(stmt, 8, seq_normal);
	sqlite3_bind_int(stmt, 9, a_count);
	sqlite3_bind_int(stmt, 10, c_count);
	sqlite3_bind_int(stmt, 11, g_count);
	sqlite3_bind_int(stmt, 12, t_count);
	sqlite3_bind_int(stmt, 13, n_count);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	sqlite3_exec(self->index_db, "COMMIT;", NULL, NULL, NULL);

	//create gzip random access index
	if(self->gzip_format){
		_pyfastx_build_gzip_index(self);
	}

	//get attributes
	_pyfastx_calc_fasta_attrs(self);
}

//load index from index file
void _pyfastx_load_index(pyfastx_Fasta *self){
	if(sqlite3_open(self->index_file, &self->index_db) != SQLITE_OK){
		PyErr_SetString(PyExc_ConnectionError, sqlite3_errmsg(self->index_db));
		return;
	}

	if(self->gzip_format){
		_pyfastx_load_gzip_index(self);
	}

	_pyfastx_calc_fasta_attrs(self);
}

/*
@param name str, sequence name
@param start int, one-based start position
@param end int, one-based end position
@param strand char, default +, - for reverse complement
*/
PyObject *get_sub_seq(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	char *name;
	int start;
	int end;
	char *strand = "+";
	static char* kwlist[] = {"name", "start", "end", "strand", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "sii|s", kwlist, &name, &start, &end, &strand)){
		return NULL;
	}

	sqlite3_stmt *stmt;
	char *sql = "SELECT * FROM seq WHERE sid=? LIMIT 1;";
	if(sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return NULL;
	}

	//byte start of sequence in file
	int byte_offset;

	//byte end of sequence in file
	int byte_len;

	//sequence length
	//int seq_len;

	//line length of fasta
	int line_len;

	//line end length
	int end_len;

	//is standard FASTA format that with the same line length
	int standard;

	//how many lines the sequence occupied
	int line_num;

	//number of bases in a incomplete line
	int tail_num;

	//offset for sub sequence start
	int offset;

	//read length
	int read_byte;

	sqlite3_bind_text(stmt, 1, name, strlen(name), NULL);
	sqlite3_step(stmt);
	byte_offset = sqlite3_column_int(stmt, 1);
	byte_len = sqlite3_column_int(stmt, 2);
	//seq_len = sqlite3_column_int(stmt, 3);
	line_len = sqlite3_column_int(stmt, 4);
	end_len = sqlite3_column_int(stmt, 5);
	standard = sqlite3_column_int(stmt, 6);
	sqlite3_finalize(stmt);

	//not standard FASTA format
	offset = byte_offset;
	read_byte = byte_len;

	if(standard){
		line_num = (end - start + 1) / (line_len - end_len);
		tail_num = (end - start + 1) % (line_len - end_len);
		offset = byte_offset + start + (start / (line_len - end_len)) * end_len - 1;
		read_byte = line_num * line_len + tail_num;
	}

	//read sequence
	char *buff = (char *)malloc(read_byte+1);

	if(self->gzip_format){
		zran_seek(self->gzip_index, offset, SEEK_SET, NULL);
		zran_read(self->gzip_index, buff, read_byte);
	} else {
		gzseek(self->gzfd, offset, SEEK_SET);
		gzread(self->gzfd, buff, read_byte);
	}

	buff[read_byte] = '\0';

	if (!strand){
		truncate_seq(buff, start, end);
	}

	return Py_BuildValue("s", buff);
}



PyObject *fasta_tp_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
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

	//initial kseqs
	obj->gzfd = gzopen(obj->file_name, "rb");
	obj->kseqs = kseq_init(obj->gzfd);

	//check input file is gzip or not
	obj->fd = fopen(obj->file_name, "rb");
	obj->gzip_format = is_gzip_format(obj->fd);

	//create index file path
	obj->index_file = (char *)malloc(strlen(file_name) + 4);
	strcpy(obj->index_file, file_name);
	strcat(obj->index_file, ".db");

	//if index file exists, connect to it
	if(obj->gzip_format){
		obj->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
	}
	
	//index database connection
	obj->index_db = NULL;

	//if build_index is True
	if(build_index){
		if(file_exists(obj->index_file)) {
			_pyfastx_load_index(obj);
		} else {
			_pyfastx_build_index(obj);
		}
	}
	
	return (PyObject *)obj;
}

void fasta_tp_dealloc(pyfastx_Fasta *self){
	kseq_destroy(self->kseqs);
	if(self->gzip_format){
		zran_free(self->gzip_index);
	}
	if(self->index_db){
		sqlite3_close(self->index_db);
	}
	gzclose(self->gzfd);
	fclose(self->fd);
	Py_TYPE(self)->tp_free(self);
}

PyObject *fasta_tp_iter(pyfastx_Fasta *self){
	kseq_rewind(self->kseqs);
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *fasta_tp_next(pyfastx_Fasta *self){
	if(kseq_read(self->kseqs) >= 0){
		if(self->uppercase){
			upper_string(self->kseqs->seq.s);
		}
		return Py_BuildValue("(ss)", self->kseqs->name.s, self->kseqs->seq.s);
	}
	return NULL;
}

PyObject *fasta_build_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	if(file_exists(self->index_file)) {
		_pyfastx_load_index(self);
	} else {
		_pyfastx_build_index(self);
	}
	return Py_BuildValue("");
}

PyObject *fasta_rebuild_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	if(file_exists(self->index_file)){
		remove(self->index_file);
	}
	_pyfastx_build_index(self);
	return Py_BuildValue("");
}

int fasta_get_item(pyfastx_Fasta *self, PyObject *key){
	if(PyUnicode_Check(key)){
		return 1;
	}
	return 0;
}

PyObject *fasta_get_key(pyfastx_Fasta *self, PyObject *key){
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

int fasta_get_len(pyfastx_Fasta *self){
	return self->seq_counts;
}

int fasta_get_val(pyfastx_Fasta *self, PyObject *key, PyObject *val){ 
	//char *name = PyUnicode_AsUTF8(key);
	return 1;
}

static PyMemberDef fasta_members[] = {
	{"file_name", T_STRING, offsetof(pyfastx_Fasta, file_name), READONLY},
	{"size", T_LONG, offsetof(pyfastx_Fasta, seq_length), READONLY},
	{"count", T_INT, offsetof(pyfastx_Fasta, seq_counts), READONLY},
	{"gc_content", T_FLOAT, offsetof(pyfastx_Fasta, gc_content), READONLY},
	{"composition", T_OBJECT, offsetof(pyfastx_Fasta, composition), READONLY},
	{NULL}
};

static PyMethodDef fasta_methods[] = {
	{"build_index", (PyCFunction)fasta_build_index, METH_VARARGS},
	{"rebuild_index", (PyCFunction)fasta_rebuild_index, METH_VARARGS},
	{"get_sub_seq", (PyCFunction)get_sub_seq, METH_VARARGS},
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

static PyMappingMethods map_methods = {
	(lenfunc)fasta_get_len,
	(binaryfunc)fasta_get_key,
	(objobjargproc)fasta_get_val,
};

PyTypeObject pyfastx_FastaType = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "Fasta",                        /* tp_name */
    sizeof(pyfastx_Fasta),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)fasta_tp_dealloc,   /* tp_dealloc */
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
    (getiterfunc)fasta_tp_iter,     /* tp_iter */
    (iternextfunc)fasta_tp_next,    /* tp_iternext */
    fasta_methods,                  /* tp_methods */
    fasta_members,                  /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    fasta_tp_new,                   /* tp_new */
};
