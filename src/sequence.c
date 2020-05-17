#include "sequence.h"
#include "structmember.h"

char *pyfastx_sequence_get_subseq(pyfastx_Sequence* self) {
	if (!self->normal || (self->parent_len == self->end && self->start == 1)) {
		pyfastx_index_get_full_seq(self->index, self->id);
	}
	
	if ((self->id == self->index->cache_chrom) && (self->start==self->index->cache_start) && (self->end==self->index->cache_end)){
		return self->index->cache_seq;
	}

	if ((self->id == self->index->cache_chrom) && (self->start>=self->index->cache_start) && (self->end<=self->index->cache_end)){
		char *buff = (char *)malloc(self->seq_len + 1);
		memcpy(buff, self->index->cache_seq + (self->start - self->index->cache_start), self->seq_len);
		buff[self->seq_len] = '\0';
		return buff;
	}
	
	if (self->index->cache_chrom) {
		free(self->index->cache_seq);
	}
	self->index->cache_seq = (char *)malloc(self->byte_len + 1);

	//Py_BEGIN_ALLOW_THREADS
	if (self->index->gzip_format) {
		zran_seek(self->index->gzip_index, self->offset, SEEK_SET, NULL);
		zran_read(self->index->gzip_index, self->index->cache_seq, self->byte_len);
	} else {
		gzseek(self->index->gzfd, self->offset, SEEK_SET);
		gzread(self->index->gzfd, self->index->cache_seq, self->byte_len);
	}

	self->index->cache_seq[self->byte_len] = '\0';

	if (self->index->uppercase) {
		remove_space_uppercase(self->index->cache_seq);
	} else {
		remove_space(self->index->cache_seq);
	}

	//Py_END_ALLOW_THREADS
	self->index->cache_chrom = self->id;
	self->index->cache_start = self->start;
	self->index->cache_end = self->end;

	return self->index->cache_seq;
}


/*PyObject *pyfastx_sequence_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	pyfastx_Sequence *obj = (pyfastx_Sequence *)type->tp_alloc(type, 0);
	return (PyObject *)obj;
}*/

void pyfastx_sequence_dealloc(pyfastx_Sequence* self) {
	free(self->name);

	if (self->desc) {
		free(self->desc);
	}

	if (self->line.l > 0) {
		free(self->line.s);
	}

	if (self->line_cache) {
		free(self->line_cache);
	}

	Py_TYPE(self)->tp_free(self);
}

PyObject *pyfastx_sequence_iter(pyfastx_Sequence* self){
	if (self->start != 1 || self->end != self->parent_len) {
		PyErr_SetString(PyExc_RuntimeError, "sliced sequence cannot be read line by line");
		return NULL;
	}

	if (self->index->gzip_format){
		zran_seek(self->index->gzip_index, self->offset, SEEK_SET, NULL);
	} else {
		gzseek(self->index->gzfd, self->offset, SEEK_SET);
	}

	if (!self->line_cache) {
		self->line_cache = (char *)malloc(1048576 + 1);
	}

	self->cache_pos = NULL;

	if (!self->line.m) {
		self->line.m = 1;
		self->line.l = 0;
		self->line.s = (char *)malloc(1);
	}

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_sequence_next(pyfastx_Sequence* self){
	uint32_t rlen, len;
	char *ret;
	
	if (self->line.l > 0) {
		self->line.l = 0;
	}

	while (1) {
		if (!self->cache_pos) {
			if (self->index->gzip_format) {
				rlen = zran_read(self->index->gzip_index, self->line_cache, 1048576);
			} else {
				rlen = gzread(self->index->gzfd, self->line_cache, 1048576);
			}

			if (rlen <= 0) {
				if (self->line.l > 0) {
					return Py_BuildValue("s", self->line.s);
				}
				return NULL;
			}

			self->line_cache[rlen] = '\0';
			self->cache_pos = self->line_cache;
		}

		ret = strchr(self->cache_pos, '\n');

		if (ret) {
			len = ret - self->cache_pos + 1;
		} else {
			len = self->line_cache + strlen(self->line_cache) - self->cache_pos;
		}

		if (len + self->line.l > self->line.m) {
			self->line.m = len + self->line.l + 1;
			self->line.s = (char *)realloc(self->line.s, self->line.m);
		}

		memcpy(self->line.s+self->line.l, self->cache_pos, len);
		self->line.l += len;
		self->line.s[self->line.l] = '\0';

		if (self->line.s[0] == '>') {
			free(self->line_cache);
			self->line_cache = NULL;
			return NULL;
		}

		if (ret) {
			self->cache_pos = ret + 1;

			if (*self->cache_pos == '\0') {
				self->cache_pos = NULL;
			}

			self->line.s[self->line.l-self->end_len] = '\0';
			
			return Py_BuildValue("s", self->line.s);
		} else {
			self->cache_pos = NULL;
		}
	}

	free(self->line_cache);
	self->line_cache = NULL;
	return NULL;
}

uint32_t pyfastx_sequence_length(pyfastx_Sequence* self){
	return self->seq_len;
}

uint16_t pyfastx_sequence_contains(pyfastx_Sequence *self, PyObject *key){
	char *seq;
	char *subseq;

	if (!PyUnicode_CheckExact(key)) {
		return 0;
	}

	seq = pyfastx_sequence_get_subseq(self);
	subseq = (char *)PyUnicode_AsUTF8(key);
	
	if(strstr(seq, subseq) != NULL){
		return 1;
	}
	
	return 0;
}

PyObject *pyfastx_sequence_get_name(pyfastx_Sequence* self, void* closure){
	if(self->start == 1 && self->end == self->parent_len){
		return Py_BuildValue("s", self->name);
	} else {
		return PyUnicode_FromFormat("%s:%d-%d", self->name, self->start, self->end);
	}
}

PyObject *pyfastx_sequence_description(pyfastx_Sequence* self, void* closure){
	sqlite3_stmt *stmt;
	int nbytes;
	int ret;
	int64_t new_offset;

	if (!self->desc) {
		const char *sql = "SELECT dlen FROM seq WHERE ID=? LIMIT 1";
		
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
			sqlite3_bind_int(stmt, 1, self->id);
			ret = sqlite3_step(stmt);
		);

		if (ret != SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
			PyErr_SetString(PyExc_RuntimeError, "can not get sequence description");
			return NULL;
		}

		PYFASTX_SQLITE_CALL(
			nbytes = sqlite3_column_int(stmt, 0);
			sqlite3_finalize(stmt);
		);

		self->desc = (char *)malloc(nbytes + 1);
		new_offset = self->offset - nbytes - self->end_len;
		
		if (self->index->gzip_format) {
			zran_seek(self->index->gzip_index, new_offset, SEEK_SET, NULL);
			zran_read(self->index->gzip_index, self->desc, nbytes);
		} else {
			gzseek(self->index->gzfd, new_offset, SEEK_SET);
			gzread(self->index->gzfd, self->desc, nbytes);
		}

		self->desc[nbytes] = '\0';
	}
	return Py_BuildValue("s", self->desc);
}

char *pyfastx_sequence_acquire(pyfastx_Sequence* self){
	char *ret;
	char *seq;

	seq = pyfastx_sequence_get_subseq(self);
	
	ret = (char *)malloc(self->seq_len + 1);
	strcpy(ret, seq);
	
	if ((self->id == self->index->cache_chrom) && (self->start>=self->index->cache_start) && (self->end<=self->index->cache_end)){
		if (!(self->start==self->index->cache_start && self->end<=self->index->cache_end)) {
			free(seq);
		}
	}
	
	return ret;
}

PyObject *pyfastx_sequence_raw(pyfastx_Sequence* self, void* closure) {
	sqlite3_stmt *stmt;
	int nbytes;
	int64_t new_offset;
	int64_t new_bytelen;
	char *buff;
	int ret;
	PyObject* retval;

	const char *sql = "SELECT dlen FROM seq WHERE ID=? LIMIT 1";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, self->id);
	ret = sqlite3_step(stmt);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(nbytes = sqlite3_column_int(stmt, 0));
	} else {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_SetString(PyExc_RuntimeError, "get sequence description length error");
		return NULL;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	new_offset = self->offset - nbytes - self->end_len - 1;
	new_bytelen = self->byte_len + nbytes + self->end_len + 1;

	if (self->parent_len == self->end && self->start == 1) {
		buff = (char *)malloc(new_bytelen + 1);

		if (self->index->gzip_format) {
			zran_seek(self->index->gzip_index, new_offset, SEEK_SET, NULL);
			zran_read(self->index->gzip_index, buff, new_bytelen);
		} else {
			gzseek(self->index->gzfd, new_offset, SEEK_SET);
			gzread(self->index->gzfd, buff, new_bytelen);
		}
		buff[new_bytelen] = '\0';
		retval = Py_BuildValue("s", buff);
	} else {
		buff = pyfastx_sequence_get_subseq(self);
		retval = PyUnicode_FromFormat(">%s:%ld-%ld\n%s\n", self->name, self->start, self->end, buff);
	}

	free(buff);
	return retval;
}

PyObject *pyfastx_sequence_seq(pyfastx_Sequence* self, void* closure){
	char *seq;
	PyObject *ret;
	seq = pyfastx_sequence_get_subseq(self);
	ret = Py_BuildValue("s", seq);

	if ((self->id == self->index->cache_chrom) && (self->start>=self->index->cache_start) && (self->end<=self->index->cache_end)){
		if (!(self->start==self->index->cache_start && self->end<=self->index->cache_end)) {
			free(seq);
		}
	}

	return ret;
}

PyObject *pyfastx_sequence_reverse(pyfastx_Sequence* self, void* closure){
	char *seq;
	PyObject *ret;

	seq = pyfastx_sequence_acquire(self);
	reverse_seq(seq);
	
	ret = Py_BuildValue("s", seq);
	free(seq);

	return ret;
}

PyObject *pyfastx_sequence_complement(pyfastx_Sequence* self, void* closure){
	char *seq;
	PyObject *ret;

	seq = pyfastx_sequence_acquire(self);
	complement_seq(seq);

	ret = Py_BuildValue("s", seq);
	free(seq);

	return ret;
}

//complement reverse sequence
PyObject *pyfastx_sequence_antisense(pyfastx_Sequence* self, void* closure){
	char *seq;
	PyObject *ret;

	seq = pyfastx_sequence_acquire(self);
	reverse_complement_seq(seq);

	ret = Py_BuildValue("s", seq);
	free(seq);

	return ret;
}

PyObject *pyfastx_sequence_repr(pyfastx_Sequence* self){
	if(self->start == 1 && self->end == self->parent_len){
		return PyUnicode_FromFormat("<Sequence> %s with length of %d", self->name, self->seq_len);
	} else {
		return PyUnicode_FromFormat("<Sequence> %s from %d to %d", self->name, self->start, self->end);
	}
}

PyObject *pyfastx_sequence_str(pyfastx_Sequence* self){
	return pyfastx_sequence_seq(self, NULL);
}

PyObject *pyfastx_seqeunce_subscript(pyfastx_Sequence* self, PyObject* item){
	if (PyIndex_Check(item)) {
		Py_ssize_t i;
		char *sub_seq;

		i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		
		if (i == -1 && PyErr_Occurred()){
			return NULL;
		}

		if (i < 0) {
			i += self->seq_len;
		}

		sub_seq = pyfastx_sequence_get_subseq(self);
		return Py_BuildValue("C", *(sub_seq+i));
	
	} else if (PySlice_Check(item)) {
		Py_ssize_t slice_start, slice_stop, slice_step, slice_len;
		pyfastx_Sequence *seq;

		if (PySlice_GetIndicesEx(item, self->seq_len, &slice_start, &slice_stop, &slice_step, &slice_len) < 0) {
			return NULL;
		}

		if (slice_len <= 0) {
			Py_RETURN_NONE;
		}

		if (slice_step == 0) {
			PyErr_SetString(PyExc_ValueError, "slice step cannot be zero");
		}

		if (slice_step != 1) {
			Py_RETURN_NONE;
		}

		//create a new sequence
		seq = (pyfastx_Sequence *)PyObject_CallObject((PyObject *)&pyfastx_SequenceType, NULL);
		if (!seq) {
			return NULL;
		}
		
		seq->start = slice_start + self->start;
		seq->end = slice_stop + self->start - 1;
		seq->id = self->id;
		seq->name = (char *)malloc(strlen(self->name) + 1);
		strcpy(seq->name, self->name);
		seq->seq_len = slice_stop - slice_start;
		seq->parent_len = self->parent_len;
		seq->line_len = self->line_len;
		seq->end_len = self->end_len;
		seq->normal = self->normal;
		seq->offset = self->offset;
		seq->byte_len = self->byte_len;
		seq->index = self->index;
		seq->line_cache = NULL;
		seq->cache_pos = NULL;
		seq->line.s = NULL;
		seq->line.l = 0;
		seq->line.m = 0;

		if (self->normal) {
			//number of the lines before slice start
			int before_sline = slice_start/(self->line_len - self->end_len);

			//number of the lines before slice stop
			int before_eline = (slice_stop + 1)/(self->line_len - self->end_len);

			//number of the lines crossed by sliced sequence
			int cross_line = before_eline - before_sline;

			//int32_t line_num = (slice_start + 1)/(self->line_len - self->end_len);
			//seq->offset = self->offset + slice_start + self->end_len*line_num;
			//seq->byte_len = seq->seq_len + seq->seq_len/self->line_len*self->end_len;
			seq->offset = self->offset + slice_start + self->end_len*before_sline;
			seq->byte_len = seq->seq_len + cross_line*self->end_len;
		}

		//Py_INCREF(seq);
		return (PyObject *)seq;
	} else {
		return NULL;
	}

}

PyObject *pyfastx_sequence_search(pyfastx_Sequence *self, PyObject *args, PyObject *kwargs){
	char* keywords[] = {"subseq", "strand", NULL};

	char *subseq;
	int sublen;
	char *seq;
	char *result;
	uint32_t start;
	int strand = '+';
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|C", keywords, &subseq, &sublen, &strand)){
		return NULL;
	}

	if (strand == '-') {
		reverse_complement_seq(subseq);
	}

	seq = pyfastx_sequence_get_subseq(self);

	result = strstr(seq, subseq);
	if(result == NULL){
		Py_RETURN_NONE;
	}
	if(strand == '-'){
		start = result - seq + sublen;
	} else {
		start = result - seq + 1;
	}
	
	return Py_BuildValue("I", start);
}

PyObject *pyfastx_sequence_gc_content(pyfastx_Sequence *self, void* closure) {
	int64_t a = 0, c = 0, g = 0, t = 0;
	int ret;
	sqlite3_stmt *stmt;

	const char *sql = "SELECT a, c, g, t FROM comp WHERE ID=? LIMIT 1";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int(stmt, 1, self->id);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW && self->start == 1 && self->end == self->seq_len) {
		PYFASTX_SQLITE_CALL(
			a = sqlite3_column_int(stmt, 0);
			c = sqlite3_column_int(stmt, 1);
			g = sqlite3_column_int(stmt, 2);
			t = sqlite3_column_int(stmt, 3);
		);
	} else {
		char *seq;
		uint32_t i;
		seq = pyfastx_sequence_get_subseq(self);

		for (i = 0; i < self->seq_len; i++) {
			switch (seq[i]) {
				case 65: case 97: ++a; break;
				case 84: case 116: ++t; break;
				case 71: case 103: ++g; break;
				case 67: case 99: ++c; break;
			}
		}
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return Py_BuildValue("f", (float)(g+c)/(a+c+g+t)*100);
}

PyObject *pyfastx_sequence_gc_skew(pyfastx_Sequence *self, void* closure) {
	int64_t c = 0, g = 0;
	int ret;
	sqlite3_stmt *stmt;

	const char *sql = "SELECT c, g FROM comp WHERE ID=? LIMIT 1";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int(stmt, 1, self->id);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW && self->start == 1 && self->end == self->seq_len) {
		PYFASTX_SQLITE_CALL(
			c = sqlite3_column_int(stmt, 0);
			g = sqlite3_column_int(stmt, 1);
		);
	} else {
		char *seq;
		uint32_t i;
		seq = pyfastx_sequence_get_subseq(self);
		for (i = 0; i < self->seq_len; i++) {
			switch (seq[i]) {
				case 71: case 103: ++g; break;
				case 67: case 99: ++c; break;
			}
		}
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return Py_BuildValue("f", (float)(g-c)/(g+c));
}

PyObject *pyfastx_sequence_composition(pyfastx_Sequence *self, void* closure) {
	sqlite3_stmt *stmt;
	int16_t i;
	int64_t c;
	const char *sql = "SELECT * FROM comp WHERE ID=?";
	PyObject *d;
	int ret;

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int(stmt, 1, self->id);
		ret = sqlite3_step(stmt);
	);

	d = PyDict_New();
	
	if (ret == SQLITE_ROW && self->start == 1 && self->end == self->seq_len) {
		for (i = 1; i < 27; i++) {
			PYFASTX_SQLITE_CALL(c = sqlite3_column_int64(stmt, i));
			if (c > 0) {
				PyDict_SetItem(d, Py_BuildValue("C", i+64), Py_BuildValue("i", c));
			}
		}
	} else {
		char *seq;
		int seq_comp[26] = {0};
		seq = pyfastx_sequence_get_subseq(self);
		for (c = 0; c < self->seq_len; c++) {
			++seq_comp[seq[c]-65];
		}

		for (i = 0; i < 26; i++) {
			c = seq_comp[i];
			if (c > 0) {
				PyDict_SetItem(d, Py_BuildValue("C", i+65), Py_BuildValue("i", c));
			}
		}
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return d;
}

static PyMappingMethods pyfastx_sequence_as_mapping = {
	(lenfunc)pyfastx_sequence_length,
	(binaryfunc)pyfastx_seqeunce_subscript,
	0
};

static PySequenceMethods pyfastx_sequence_as_sequence = {
	0, /*sq_length*/
	0, /*sq_concat*/
	0, /*sq_repeat*/
	0, /*sq_item*/
	0, /*sq_slice */
	0, /*sq_ass_item*/
	0, /*sq_ass_splice*/
	(objobjproc)pyfastx_sequence_contains, /*sq_contains*/
	0, /*sq_inplace_concat*/
	0, /*sq_inplace_repeat*/
};

static PyMethodDef pyfastx_sequence_methods[] = {
	{"search", (PyCFunction)pyfastx_sequence_search, METH_VARARGS|METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyGetSetDef pyfastx_sequence_getsets[] = {
	{"name", (getter)pyfastx_sequence_get_name, NULL, NULL, NULL},
	{"raw", (getter)pyfastx_sequence_raw, NULL, NULL, NULL},
	{"seq", (getter)pyfastx_sequence_seq, NULL, NULL, NULL},
	{"reverse", (getter)pyfastx_sequence_reverse, NULL, NULL, NULL},
	{"complement", (getter)pyfastx_sequence_complement, NULL, NULL, NULL},
	{"antisense", (getter)pyfastx_sequence_antisense, NULL, NULL, NULL},
	{"description", (getter)pyfastx_sequence_description, NULL, NULL, NULL},
	{"gc_content", (getter)(pyfastx_sequence_gc_content), NULL, NULL, NULL},
	{"gc_skew", (getter)(pyfastx_sequence_gc_skew), NULL, NULL, NULL},
	{"composition", (getter)(pyfastx_sequence_composition), NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef pyfastx_sequence_members[] = {
	//{"name", T_STRING, offsetof(pyfastx_Sequence, name), READONLY},
	{"id", T_INT, offsetof(pyfastx_Sequence, id), READONLY},
	{"start", T_INT, offsetof(pyfastx_Sequence, start), READONLY},
	{"end", T_INT, offsetof(pyfastx_Sequence, end), READONLY},
	//{"length", T_INT, offsetof(pyfastx_Sequence, seq_len), READONLY},
	{NULL}
};

PyTypeObject pyfastx_SequenceType = {
    //PyVarObject_HEAD_INIT(&PyType_Type, 0)
    PyVarObject_HEAD_INIT(NULL, 0)
    "Sequence",                     /* tp_name */
    sizeof(pyfastx_Sequence),       /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_sequence_dealloc,                              /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)pyfastx_sequence_repr,                              /* tp_repr */
    0,                              /* tp_as_number */
    &pyfastx_sequence_as_sequence,                              /* tp_as_sequence */
    &pyfastx_sequence_as_mapping,   /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    (reprfunc)pyfastx_sequence_str,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,             /* tp_flags */
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    (getiterfunc)pyfastx_sequence_iter,                              /* tp_iter */
    (iternextfunc)pyfastx_sequence_next,                              /* tp_iternext */
    pyfastx_sequence_methods,       /* tp_methods */
    pyfastx_sequence_members,       /* tp_members */
    pyfastx_sequence_getsets,       /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    PyType_GenericNew,           /* tp_new */
};
