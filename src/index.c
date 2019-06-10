#include "fasta.h"
#include "util.h"
//#include "structmember.h"
KSEQ_INIT(gzFile, gzread, gzrewind)
/*
static int use_index(FastxObject *self){
	return 0;
}*/

static void build_gzip_index(FastxObject *self){
	sqlite3_stmt *stmt;
	zran_init(self->gzip_index, self->fd, 0, 0, 0, ZRAN_AUTO_BUILD);
	zran_build_index(self->gzip_index, 0, 0);

	//create temp gzip index file
	char *temp_index = malloc(strlen(self->file_name) + 5);
	strcpy(temp_index, self->file_name);
	strcat(temp_index, ".tmp");
	FILE* fd = fopen(temp_index, "wb");
	zran_export_index(self->gzip_index, fd);
	
	long fsize = ftell(fd);
	fseek(fd, 0, SEEK_SET);
	char *buff = malloc(fsize + 1);
	int ret = fread(buff, 1, fsize, fd);
	fclose(fd);
	remove(temp_index);

	sqlite3_prepare_v2(self->db, "INSERT INTO gzindex VALUES (NULL, ?)", -1, &stmt, NULL);
	sqlite3_bind_blob(stmt, 1, buff, strlen(buff), NULL);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
}

static void load_gzip_index(FastxObject *self){
	sqlite3_stmt *stmt;
	zran_init(self->gzip_index, self->fd, 0, 0, 0, ZRAN_AUTO_BUILD);
	char *temp_index = malloc(strlen(self->file_name) + 5);
	strcpy(temp_index, self->file_name);
	strcat(temp_index, ".tmp");
	FILE* fd = fopen(temp_index, "wb+");

	sqlite3_prepare_v2(self->db, "SELECT content FROM gzindex;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	const char *buff = sqlite3_column_blob(stmt, 0);
	fwrite(buff, 1, strlen(buff), fd);
	fclose(fd);
	remove(temp_index);
}


PyObject *build_index(FastxObject *self, PyObject *args, PyObject *kwargs){
	// 1 force rebuild index, 0 for not rebuild index
	int force = 0;
	
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

	//kwargs parameters
	static char* kwlist[] = {"force", NULL};

	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "|p", kwlist, &force)){
		return NULL;
	}

	//check index file exists, if exists do not build index
	if(file_exists(self->index_file)){
		if(force){
			remove(self->index_file);
		} else {
			return Py_BuildValue("i", 1);
		}
	}

	ret = sqlite3_open(self->index_file, &self->db);
	if(ret != SQLITE_OK){
		return Py_BuildValue("i", 0);
	}

	//create index database
	const char *create_sql = "\
		CREATE TABLE seq (\
			sid TEXTPRIMARY KEY,\
			offset INTEGER,\
			blen INTEGER,\
			slen INTEGER,\
			llen INTEGER,\
			elen INTEGER,\
			norm INTEGER,\
			a INTEGER,\
			c INTEGER,\
			g INTEGER,\
			t INTEGER,\
			n INTEGER\
		);\
		CREATE TABLE gzindex (\
			ID INTEGER PRIMARY KEY,\
			content BLOB\
		);";

	ret = sqlite3_exec(self->db, create_sql, NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "create table error");
	}

	ret = sqlite3_exec(self->db, "PRAGMA synchronous=OFF", NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "pragma");
	}

	ret = sqlite3_exec(self->db, "begin", NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "begin error");
	}

	char *insert_sql = "INSERT INTO seq VALUES (?,?,?,?,?,?,?,?,?,?,?,?)";
	
	sqlite3_prepare_v2(self->db, insert_sql, -1, &stmt, NULL);

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
				sqlite3_bind_text(stmt, 1, self->kseqs->name.s, self->kseqs->name.l, NULL);
				sqlite3_bind_int(stmt, 2, start);
				sqlite3_bind_int(stmt, 3, position-start-1);
				sqlite3_bind_int(stmt, 4, seq_len);
				sqlite3_bind_int(stmt, 5, line_len);
				sqlite3_bind_int(stmt, 6, line_end);
				sqlite3_bind_int(stmt, 7, seq_normal);
				sqlite3_bind_int(stmt, 8, a_count);
				sqlite3_bind_int(stmt, 9, c_count);
				sqlite3_bind_int(stmt, 10, g_count);
				sqlite3_bind_int(stmt, 11, t_count);
				sqlite3_bind_int(stmt, 12, n_count);
				sqlite3_step(stmt);
				sqlite3_reset(stmt);
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

	sqlite3_bind_text(stmt, 1, self->kseqs->name.s, self->kseqs->name.l, NULL);
	sqlite3_bind_int(stmt, 2, start);
	sqlite3_bind_int(stmt, 3, position-start-1);
	sqlite3_bind_int(stmt, 4, seq_len);
	sqlite3_bind_int(stmt, 5, line_len);
	sqlite3_bind_int(stmt, 6, line_end);
	sqlite3_bind_int(stmt, 7, seq_normal);
	sqlite3_bind_int(stmt, 8, a_count);
	sqlite3_bind_int(stmt, 9, c_count);
	sqlite3_bind_int(stmt, 10, g_count);
	sqlite3_bind_int(stmt, 11, t_count);
	sqlite3_bind_int(stmt, 12, n_count);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	sqlite3_exec(self->db, "commit;", NULL, NULL, NULL);

	//create gzip random access index
	if(self->gzip){
		/*
		zran_init(self->gzip_index, self->fd, 0, 0, 0, ZRAN_AUTO_BUILD);
		zran_build_index(self->gzip_index, 0, 0);

		//create temp gzip index file
		char *temp_index = malloc(strlen(self->file_name) + 5);
		strcpy(temp_index, self->file_name);
		strcat(temp_index, ".tmp");
		FILE* fd = fopen(temp_index, "wb");
		zran_export_index(self->gzip_index, fd);
		
		long fsize = ftell(fd);
		fseek(fd, 0, SEEK_SET);
		char *buff = malloc(fsize + 1);
		ret = fread(buff, 1, fsize, fd);
		fclose(fd);
		remove(temp_index);

		sqlite3_prepare_v2(self->db, "INSERT INTO gzindex VALUES (NULL, ?)", -1, &stmt, NULL);
		sqlite3_bind_blob(stmt, 1, buff, strlen(buff), NULL);
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);*/
		build_gzip_index(self);
	}
	
	return Py_BuildValue("i", 1);
}

PyObject* test(FastxObject *self, PyObject *args){
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
}

/*
@param name str, sequence name
@param start int, one-based start position
@param end int, one-based end position
@param strand char, default +, - for reverse complement
*/
PyObject *get_sub_seq(FastxObject *self, PyObject *args, PyObject *kwargs){
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
	if(sqlite3_prepare_v2(self->db, sql, -1, &stmt, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->db));
		return NULL;
	}

	//byte start of sequence in file
	int byte_offset;

	//byte end of sequence in file
	int byte_len;

	//sequence length
	int seq_len;

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
	seq_len = sqlite3_column_int(stmt, 3);
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

	if(self->gzip){
		zran_seek(self->gzip_index, offset, SEEK_SET, NULL);
		zran_read(self->gzip_index, buff, read_byte);
	} else {
		gzseek(self->gzfp, offset, SEEK_SET);
		gzread(self->gzfp, buff, read_byte);
	}

	buff[read_byte] = '\0';

	if (!strand){
		truncate_seq(buff, start, end);
	}

	return Py_BuildValue("s", buff);
}

PyObject *fastx_tp_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	char *file_name;
	int capital = 0;
	
	static char* keywords[] = {"file_name", "capital", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s|p", keywords, &file_name, &capital)){
		return NULL;
	}

	//check input sequence file is whether exists
	if(!file_exists(file_name)){
		return PyErr_Format(PyExc_FileExistsError, "input sequence file %s does not exists", file_name);
	}

	//create Fastx class
	FastxObject *obj = (FastxObject *)type->tp_alloc(type, 0);
	if (!obj){
		return NULL;
	}
	
	//initial sequence file name
	obj->file_name = (char *)malloc(strlen(file_name)+1);
	strcpy(obj->file_name, file_name);

	obj->capital = capital;

	//initial kseqs
	obj->gzfp = gzopen(obj->file_name, "rb");
	obj->kseqs = kseq_init(obj->gzfp);

	//open file and check file is gzip
	obj->fd = fopen(obj->file_name, "rb");
	obj->gzip = is_gzip(obj->fd);

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
		if(self->capital){
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
