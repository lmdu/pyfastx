#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "sequence.h"
#include "structmember.h"

void pyfastx_sequence_continue_read(pyfastx_Sequence* self) {
	int header_len;

	Py_ssize_t offset;
	Py_ssize_t bytelen;

	//current file read location
	Py_ssize_t current;
	Py_ssize_t gap;
	Py_ssize_t rlen;

	//raw string offset and bytelen
	header_len = self->desc_len + self->end_len + 1;
	offset = self->offset - header_len;
	bytelen = self->byte_len + header_len;

	self->raw = (char *)malloc(bytelen + 1);

	current = gztell(self->index->gzfd);
	gap = offset - current;

	if (self->index->gzip_format) {
		if (gap >= 0) {
			while (gap > 0) {
				rlen = gap > bytelen ? bytelen : gap;
				gzread(self->index->gzfd, self->raw, rlen);
				gap -= rlen;
			}
			gzread(self->index->gzfd, self->raw, bytelen);
		} else {
			zran_seek(self->index->gzip_index, offset, SEEK_SET, NULL);
			zran_read(self->index->gzip_index, self->raw, bytelen);
		}
	} else {
		if (gap != 0) {
			gzseek(self->index->gzfd, offset, SEEK_SET);
		}
		gzread(self->index->gzfd, self->raw, bytelen);
	}
	self->raw[bytelen] = '\0';

	/*if (gap > 0 && self->index->gzip_format) {
		while (gap > 0) {
			rlen = gap > bytelen ? bytelen : gap;
			gzread(self->index->gzfd, self->raw, rlen);
			gap -= rlen;
		}
	} else {
		gzseek(self->index->gzfd, offset, SEEK_SET);
	}

	gzread(self->index->gzfd, self->raw, bytelen);
	self->raw[bytelen] = '\0';*/

	self->desc = (char *)malloc(self->desc_len + 1);
	memcpy(self->desc, self->raw+1, self->desc_len);
	self->desc[self->desc_len] = '\0';

	//copy sequence to cache
	if (self->byte_len >= self->index->cache_seq.m) {
		self->index->cache_seq.m = self->byte_len + 1;
		self->index->cache_seq.s = (char *)realloc(self->index->cache_seq.s, self->index->cache_seq.m);
	}

	memcpy(self->index->cache_seq.s, self->raw+self->desc_len+self->end_len+1, self->byte_len);
	self->index->cache_seq.s[self->byte_len] = '\0';

	if (self->index->uppercase) {
		self->index->cache_seq.l = remove_space_uppercase(self->index->cache_seq.s, self->byte_len);
	} else {
		self->index->cache_seq.l = remove_space(self->index->cache_seq.s, self->byte_len);
	}

	self->index->cache_chrom = self->id;
	self->index->cache_start = 1;
	self->index->cache_end = self->seq_len;
	self->index->cache_full = 1;
}

char *pyfastx_sequence_get_fullseq(pyfastx_Sequence* self) {
	if ((self->id == self->index->cache_chrom) && self->index->cache_full) {
		return self->index->cache_seq.s;
	}

	if (strlen(self->name) >= self->index->cache_name.m) {
		self->index->cache_name.m = strlen(self->name) + 1;
		self->index->cache_name.s = (char *)realloc(self->index->cache_name.s, self->index->cache_name.m);
	}

	strcpy(self->index->cache_name.s, self->name);

	pyfastx_index_fill_cache(self->index, self->offset, self->byte_len);

	self->index->cache_chrom = self->id;
	self->index->cache_start = 1;
	self->index->cache_end = self->seq_len;
	self->index->cache_full = 1;

	return self->index->cache_seq.s;
}

char *pyfastx_sequence_get_subseq(pyfastx_Sequence* self) {
	if (self->complete || !self->normal) {
		pyfastx_sequence_get_fullseq(self);
	}
	
	if ((self->id == self->index->cache_chrom) && (self->start==self->index->cache_start) && (self->end==self->index->cache_end)){
		return self->index->cache_seq.s;
	}

	if ((self->id == self->index->cache_chrom) && (self->start>=self->index->cache_start) && (self->end<=self->index->cache_end)){
		return self->index->cache_seq.s + (self->start - self->index->cache_start);
	}

	if (self->index->cache_name.s) {
		self->index->cache_name.s[0] = '\0';
	}

	pyfastx_index_fill_cache(self->index, self->offset, self->byte_len);

	//Py_END_ALLOW_THREADS
	self->index->cache_chrom = self->id;
	self->index->cache_start = self->start;
	self->index->cache_end = self->end;
	self->index->cache_full = 0;

	return self->index->cache_seq.s;
}

void pyfastx_sequence_dealloc(pyfastx_Sequence* self) {
	free(self->name);

	if (self->desc) {
		free(self->desc);
	}

	if (self->raw) {
		free(self->raw);
	}

	if (self->line.l > 0) {
		free(self->line.s);
	}

	if (self->line_cache) {
		free(self->line_cache);
	}

	Py_DECREF(self->index->fasta);

	self->index = NULL;
	self->cache_pos = NULL;

	Py_TYPE(self)->tp_free((PyObject *)self);
}

void pyfastx_sequence_free_subseq(pyfastx_Sequence* self, char *seq) {
	if ((self->id == self->index->cache_chrom) && (self->start>=self->index->cache_start) && (self->end<=self->index->cache_end)){
		if (!(self->start==self->index->cache_start && self->end==self->index->cache_end)) {
			free(seq);
		}
	}
}

PyObject *pyfastx_sequence_iter(pyfastx_Sequence* self){
	if (!self->complete) {
		PyErr_SetString(PyExc_RuntimeError, "sliced subsequence cannot be read line by line");
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

PyObject *pyfastx_sequence_next(pyfastx_Sequence* self) {
	char *ret;

	//read length each time
	Py_ssize_t rlen;

	//real length
	Py_ssize_t len;

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

Py_ssize_t pyfastx_sequence_length(pyfastx_Sequence* self){
	return self->seq_len;
}

int pyfastx_sequence_contains(pyfastx_Sequence *self, PyObject *key){
	char *seq;
	char *subseq;
	char *ret;

	Py_ssize_t sublen;

	if (!PyUnicode_CheckExact(key)) {
		return 0;
	}

	if (self->index->iterating) {
		pyfastx_sequence_continue_read(self);
	}

	seq = pyfastx_sequence_get_subseq(self);
	subseq = (char *)PyUnicode_AsUTF8AndSize(key, &sublen);
	ret = str_n_str(seq, subseq, sublen, self->seq_len);

	return ret ? 1 : 0;
}

PyObject *pyfastx_sequence_get_name(pyfastx_Sequence* self, void* closure){
	if(self->complete){
		return Py_BuildValue("s", self->name);
	} else {
		return PyUnicode_FromFormat("%s:%ld-%ld", self->name, self->start, self->end);
	}
}

PyObject *pyfastx_sequence_description(pyfastx_Sequence* self, void* closure){
	Py_ssize_t new_offset;

	if (self->index->iterating) {
		pyfastx_sequence_continue_read(self);
	}

	if (!self->desc) {
		self->desc = (char *)malloc(self->desc_len + 1);
		new_offset = self->offset - self->desc_len - self->end_len;
		pyfastx_index_random_read(self->index, self->desc, new_offset, self->desc_len);
	}
	return Py_BuildValue("s", self->desc);
}

PyObject *pyfastx_sequence_raw(pyfastx_Sequence* self, void* closure) {
	Py_ssize_t new_offset;
	Py_ssize_t new_bytelen;

	if (self->index->iterating) {
		pyfastx_sequence_continue_read(self);
	}

	if (!self->raw) {
		if (self->complete) {
			new_offset = self->offset - self->desc_len - self->end_len - 1;
			new_bytelen = self->byte_len + self->desc_len + self->end_len + 1;
		} else {
			new_offset = self->offset;
			new_bytelen = self->byte_len;
		}

		self->raw = (char *)malloc(new_bytelen + 1);
		pyfastx_index_random_read(self->index, self->raw, new_offset, new_bytelen);
	}
	return Py_BuildValue("s", self->raw);
}

PyObject *pyfastx_sequence_seq(pyfastx_Sequence* self, void* closure){
	char *seq;
	PyObject *ret;

	if (self->index->iterating) {
		pyfastx_sequence_continue_read(self);
	}

	seq = pyfastx_sequence_get_subseq(self);
	ret = PyUnicode_New(self->seq_len, 127);
	memcpy(PyUnicode_1BYTE_DATA(ret), seq, self->seq_len);

	return ret;
}

PyObject *pyfastx_sequence_reverse(pyfastx_Sequence* self, void* closure){
	char *seq;
	char *data;
	PyObject *ret;

	seq = pyfastx_sequence_get_subseq(self);

	ret = PyUnicode_New(self->seq_len, 127);
	data = (char *)PyUnicode_1BYTE_DATA(ret);
	memcpy(data, seq, self->seq_len);

	reverse_seq(data);

	return ret;
}

PyObject *pyfastx_sequence_complement(pyfastx_Sequence* self, void* closure){
	char *seq;
	char *data;
	PyObject *ret;

	seq = pyfastx_sequence_get_subseq(self);

	ret = PyUnicode_New(self->seq_len, 127);
	data = (char *)PyUnicode_1BYTE_DATA(ret);
	memcpy(data, seq, self->seq_len);

	complement_seq(data);

	return ret;
}

//complement reverse sequence
PyObject *pyfastx_sequence_antisense(pyfastx_Sequence* self, void* closure){
	char *seq;
	char *data;
	PyObject *ret;

	seq = pyfastx_sequence_get_subseq(self);
	ret = PyUnicode_New(self->seq_len, 127);
	data = (char *)PyUnicode_1BYTE_DATA(ret);
	memcpy(data, seq, self->seq_len);

	reverse_complement_seq(data);

	return ret;
}

PyObject *pyfastx_sequence_repr(pyfastx_Sequence* self){
	if(self->complete){
		return PyUnicode_FromFormat("<Sequence> %s with length of %ld", self->name, self->seq_len);
	} else {
		return PyUnicode_FromFormat("<Sequence> %s from %ld to %ld", self->name, self->start, self->end);
	}
}

PyObject *pyfastx_sequence_str(pyfastx_Sequence* self){
	return pyfastx_sequence_seq(self, NULL);
}

PyObject *pyfastx_sequence_subscript(pyfastx_Sequence* self, PyObject* item){
	char *sub_seq;

	int before_sline;
	int before_eline;
	int cross_line;

	Py_ssize_t i;
	Py_ssize_t slice_start;
	Py_ssize_t slice_stop;
	Py_ssize_t slice_step;
	Py_ssize_t slice_len;

	pyfastx_Sequence *seq;

	if (PyIndex_Check(item)) {
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
		if (PySlice_Unpack(item, &slice_start, &slice_stop, &slice_step) < 0) {
			return NULL;
		}

		slice_len = PySlice_AdjustIndices(self->seq_len, &slice_start, &slice_stop, slice_step);

		//if (PySlice_GetIndicesEx(item, self->seq_len, &slice_start, &slice_stop, &slice_step, &slice_len) < 0) {
		//	return NULL;
		//}

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
		//seq->parent_len = self->parent_len;
		seq->line_len = self->line_len;
		seq->end_len = self->end_len;
		seq->normal = self->normal;
		seq->offset = self->offset;
		seq->byte_len = self->byte_len;
		seq->index = self->index;
		seq->line_cache = NULL;
		seq->cache_pos = NULL;
		kstring_init(seq->line);

		Py_INCREF(self->index->fasta);

		//check sequence is complete or not
		if (self->complete && seq->seq_len == self->seq_len) {
			seq->complete = 1;
		} else {
			seq->complete = 0;
		}

		if (self->normal) {
			//number of the lines before slice start
			before_sline = slice_start/(self->line_len - self->end_len);

			//number of the lines before slice stop
			before_eline = slice_stop/(self->line_len - self->end_len);

			//number of the lines crossed by sliced sequence
			cross_line = before_eline - before_sline;

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
	int strand = '+';

	char *subseq;
	char *seq;
	char *result;

	Py_ssize_t start;
	Py_ssize_t sublen;

	char* keywords[] = {"subseq", "strand", NULL};
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|C", keywords, &subseq, &sublen, &strand)){
		return NULL;
	}

	if (strand == '-') {
		reverse_complement_seq(subseq);
	}

	if (self->index->iterating) {
		pyfastx_sequence_continue_read(self);
	}

	seq = pyfastx_sequence_get_subseq(self);

	result = str_n_str(seq, subseq, sublen, self->seq_len);

	if (result != NULL) {
		if (strand == '-') {
			start = result - seq + sublen;
		} else {
			start = result - seq + 1;
		}
	}
	
	if (result == NULL) {
		Py_RETURN_NONE;
	}
	
	return Py_BuildValue("n", start);
}

PyObject *pyfastx_sequence_gc_content(pyfastx_Sequence *self, void* closure) {
	int l;
	int ret;
	char *seq;

	sqlite3_stmt *stmt;

	Py_ssize_t i, n;
	Py_ssize_t a = 0, c = 0, g = 0, t = 0;

	const char *sql = "SELECT * FROM comp WHERE seqid=?";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int64(stmt, 1, self->id);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW && self->start == 1 && self->end == self->seq_len) {
		while (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(
				l = sqlite3_column_int(stmt, 2);
				n = sqlite3_column_int64(stmt, 3);
				ret = sqlite3_step(stmt);
			);

			switch (l) {
				case 65:
				case 97:
					a += n;
					break;

				case 67:
				case 99:
					c += n;
					break;

				case 71:
				case 103:
					g += n;
					break;

				case 84:
				case 116:
					t += n;
					break;
			}
		}
	} else {
		seq = pyfastx_sequence_get_subseq(self);

		for (i = 0; i < self->seq_len; ++i) {
			switch (seq[i]) {
				case 65:
				case 97:
					++a;
					break;

				case 84:
				case 116:
					++t;
					break;

				case 71:
				case 103:
					++g;
					break;

				case 67:
				case 99:
					++c;
					break;
			}
		}
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return Py_BuildValue("f", (float)(g+c)/(a+c+g+t)*100);
}

PyObject *pyfastx_sequence_gc_skew(pyfastx_Sequence *self, void* closure) {
	int l;
	int ret;
	char *seq;
	sqlite3_stmt *stmt;

	Py_ssize_t i, n;
	Py_ssize_t c = 0, g = 0;

	const char *sql = "SELECT * FROM comp WHERE seqid=?";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int64(stmt, 1, self->id);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW && self->start == 1 && self->end == self->seq_len) {
		while (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(
				l = sqlite3_column_int(stmt, 2);
				n = sqlite3_column_int64(stmt, 3);
				ret = sqlite3_step(stmt);
			);

			switch (l) {
				case 67:
				case 99:
					c += n;
					break;

				case 71:
				case 103:
					g += n;
					break;
			}
		}
	} else {
		seq = pyfastx_sequence_get_subseq(self);

		for (i = 0; i < self->seq_len; ++i) {
			switch (seq[i]) {
				case 67:
				case 99:
					++c;
					break;

				case 71:
				case 103:
					++g;
					break;	
			}
		}
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return Py_BuildValue("f", (float)(g-c)/(g+c));
}

PyObject *pyfastx_sequence_composition(pyfastx_Sequence *self, void* closure) {
	int i;
	int l;
	int ret;
	char *seq;

	sqlite3_stmt *stmt;
	
	Py_ssize_t n;
	Py_ssize_t seq_comp[128] = {0};

	PyObject *d;
	PyObject *b;
	PyObject *c;

	const char *sql = "SELECT * FROM comp WHERE ID=?";

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int64(stmt, 1, self->id);
		ret = sqlite3_step(stmt);
	);

	d = PyDict_New();
	
	if (ret == SQLITE_ROW && self->start == 1 && self->end == self->seq_len) {
		for (i = 1; i < 128; ++i) {
			PYFASTX_SQLITE_CALL(
				l = sqlite3_column_int(stmt, 2);
				n = sqlite3_column_int64(stmt, 3);
				ret = sqlite3_step(stmt);
			);

			if (n > 0 && l != 13) {
				b = Py_BuildValue("C", l);
				c = Py_BuildValue("n", n);
				PyDict_SetItem(d, b, c);
				Py_DECREF(b);
				Py_DECREF(c);
			}
		}
	} else {
		seq = pyfastx_sequence_get_subseq(self);

		for (i = 0; i < self->seq_len; ++i) {
			++seq_comp[(unsigned char)seq[i]];
		}

		for (l = 0; l < 128; ++l) {
			n = seq_comp[l];

			if (n > 0 && l != 13) {
				b = Py_BuildValue("C", l);
				c = Py_BuildValue("i", n);
				PyDict_SetItem(d, b, c);
				Py_DECREF(b);
				Py_DECREF(c);
			}
		}
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return d;
}

static PyMappingMethods pyfastx_sequence_as_mapping = {
	.mp_length = (lenfunc)pyfastx_sequence_length,
	.mp_subscript = (binaryfunc)pyfastx_sequence_subscript,
};

static PySequenceMethods pyfastx_sequence_as_sequence = {
	.sq_contains = (objobjproc)pyfastx_sequence_contains,
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
	{"id", T_PYSSIZET, offsetof(pyfastx_Sequence, id), READONLY},
	{"start", T_PYSSIZET, offsetof(pyfastx_Sequence, start), READONLY},
	{"end", T_PYSSIZET, offsetof(pyfastx_Sequence, end), READONLY},
	{NULL}
};

PyTypeObject pyfastx_SequenceType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Sequence",
    .tp_basicsize = sizeof(pyfastx_Sequence),
    .tp_dealloc = (destructor)pyfastx_sequence_dealloc,
    .tp_repr = (reprfunc)pyfastx_sequence_repr,
    .tp_as_sequence = &pyfastx_sequence_as_sequence,
    .tp_as_mapping = &pyfastx_sequence_as_mapping,
    .tp_str = (reprfunc)pyfastx_sequence_str,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_iter = (getiterfunc)pyfastx_sequence_iter,
    .tp_iternext = (iternextfunc)pyfastx_sequence_next,
    .tp_methods = pyfastx_sequence_methods,
    .tp_members = pyfastx_sequence_members,
    .tp_getset = pyfastx_sequence_getsets,
    .tp_alloc = PyType_GenericAlloc,
    .tp_new = PyType_GenericNew,
};
