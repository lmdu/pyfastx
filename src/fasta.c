#include "fasta.h"
#include "util.h"
#include "identifier.h"
#include "structmember.h"

/*calculate fasta attributes including sequence count, length,
composition (ATGCN count) and GC content
*/
void pyfastx_calc_fasta_attrs(pyfastx_Fasta *self){
	//ACGTN nucleotide counts
	//int64_t a, c, g, t, n;

	sqlite3_stmt *stmt;
	
	//sequence count
	sqlite3_prepare_v2(self->index->index_db, "SELECT COUNT(*) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	self->seq_counts = sqlite3_column_int(stmt, 0);
	sqlite3_reset(stmt);

	//sequence length
	sqlite3_prepare_v2(self->index->index_db, "SELECT SUM(slen) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	self->seq_length = sqlite3_column_int64(stmt, 0);
	//sqlite3_reset(stmt);

	//calculate base counts
	/*sqlite3_prepare_v2(self->index->index_db, "SELECT SUM(a),SUM(c),SUM(g),SUM(t),SUM(n) FROM seq LIMIT 1;", -1, &stmt, NULL);
	sqlite3_step(stmt);
	a = sqlite3_column_int64(stmt, 0);
	c = sqlite3_column_int64(stmt, 1);
	g = sqlite3_column_int64(stmt, 2);
	t = sqlite3_column_int64(stmt, 3);
	n = sqlite3_column_int64(stmt, 4);
	self->composition = Py_BuildValue("{s:K,s:K,s:K,s:K,s:K}", "A", a, "C", c, "G", g, "T", t, "N", n);
	
	//calc GC content
	self->gc_content = (float)(g+c)/(a+c+g+t)*100;

	//calc GC skew
	self->gc_skew = (float)(g-c)/(g+c);*/

	sqlite3_finalize(stmt);
}

PyObject *pyfastx_fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	//fasta file path
	char *file_name;

	//bool value for uppercase sequence
	int uppercase = 1;

	//build index or not
	int build_index = 1;

	//key function for seperating name
	PyObject *key_func = Py_None;

	//paramters for fasta object construction
	static char* keywords[] = {"file_name", "uppercase", "build_index", "key_func", NULL};
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s|iiO", keywords, &file_name, &uppercase, &build_index, &key_func)){
		return NULL;
	}

	if ((key_func != Py_None) && !PyCallable_Check(key_func)) {
		PyErr_SetString(PyExc_TypeError, "key_func must be a callable function");
		return NULL;
	}

	//check input sequence file is whether exists
	if(!file_exists(file_name)){
		return PyErr_Format(PyExc_FileExistsError, "input fasta file %s does not exists", file_name);
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

	//create index
	obj->index = pyfastx_init_index(obj->file_name, uppercase, key_func);

	//if build_index is True
	if (build_index) {
		pyfastx_build_index(obj->index);
		pyfastx_calc_fasta_attrs(obj);
	}
	
	return (PyObject *)obj;
}

void pyfastx_fasta_dealloc(pyfastx_Fasta *self){
	pyfastx_index_free(self->index);
	Py_TYPE(self)->tp_free(self);
}

PyObject *pyfastx_fasta_iter(pyfastx_Fasta *self){
	pyfastx_rewind_index(self->index);
	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fasta_repr(pyfastx_Fasta *self){
	return PyUnicode_FromFormat("<Fasta> %s contains %d sequences", self->file_name, self->seq_counts);
}

PyObject *pyfastx_fasta_next(pyfastx_Fasta *self){
	return pyfastx_get_next_seq(self->index);
}

PyObject *pyfastx_fasta_build_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	if(!file_exists(self->index->index_file)){
		pyfastx_build_index(self->index);
		pyfastx_calc_fasta_attrs(self);
	}
	Py_RETURN_NONE;
}

PyObject *pyfastx_fasta_rebuild_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	sqlite3_close(self->index->index_db);

	if (file_exists(self->index->index_file)) {
		remove(self->index->index_file);
	}
	pyfastx_build_index(self->index);
	pyfastx_calc_fasta_attrs(self);
	Py_RETURN_NONE;
}

PyObject *pyfastx_fasta_fetch(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	static char* keywords[] = {"name", "intervals", "strand", NULL};

	char *name;
	char *seq;
	PyObject *intervals;
	uint64_t start;
	uint64_t end;
	int strand = '+';
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "sO|C", keywords, &name, &intervals, &strand)){
		return NULL;
	}

	if(!PyTuple_Check(intervals) && !PyList_Check(intervals)){
		PyErr_SetString(PyExc_ValueError, "intervals must be list or tuple");
		return NULL;
	}

	if(PyList_Check(intervals)){
		intervals = PyList_AsTuple(intervals);
	}

	PyObject *item;

	item = PyTuple_GetItem(intervals, 0);
	Py_ssize_t size = PyTuple_Size(intervals);

	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	
	//select sql statement, chrom indicates seq name or chromomsome
	const char* sql = "SELECT ID FROM seq WHERE chrom=? LIMIT 1;";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, name, -1, NULL);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		return PyErr_Format(PyExc_NameError, "Sequence %s does not exists", name);
	}

	uint32_t chrom = sqlite3_column_int(stmt, 0);
	
	uint32_t seq_len;
	char *sub_seq;

	seq = pyfastx_index_get_full_seq(self->index, chrom);

	if(size == 2 && PyLong_Check(item)){
		start = PyLong_AsLong(item);
		end = PyLong_AsLong(PyTuple_GetItem(intervals, 1));

		if(start > end){
			PyErr_SetString(PyExc_ValueError, "start position > end position");
			return NULL;
		}

		seq_len = end - start + 1;

		sub_seq = (char *)malloc(seq_len + 1);
		memcpy(sub_seq, seq+start-1, seq_len);
		sub_seq[seq_len] = '\0';
	} else {
		uint32_t i;
		uint32_t j = 0;
		sub_seq = (char *)malloc(strlen(seq) + 1);

		for(i=0; i<size; i++){
			item = PyTuple_GetItem(intervals, i);
			if(PyList_Check(item)){
				item = PyList_AsTuple(item);
			}
			start = PyLong_AsLong(PyTuple_GetItem(item, 0));
			end = PyLong_AsLong(PyTuple_GetItem(item, 1));
			seq_len = end - start + 1;

			if(start > end){
				PyErr_SetString(PyExc_ValueError, "start position > end position");
				return NULL;
			}

			memcpy(sub_seq+j, seq+start-1, seq_len);
			j += seq_len;
		}
		sub_seq[j] = '\0';
	}

	if(strand == '-'){
		reverse_seq(sub_seq);
		complement_seq(sub_seq);
	}

	return Py_BuildValue("s", sub_seq);
}

PyObject *pyfastx_fasta_keys(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	pyfastx_Identifier *ids = PyObject_New(pyfastx_Identifier, &pyfastx_IdentifierType);
	if(!ids){
		return NULL;
	}

	ids->index_db = self->index->index_db;
	ids->stmt = NULL;
	ids->seq_counts = self->seq_counts;
	ids->sort = 1;
	ids->order = 0;

	Py_INCREF(ids);
	return (PyObject *)ids;
}

PyObject *pyfastx_fasta_subscript(pyfastx_Fasta *self, PyObject *item){
	
	if (PyIndex_Check(item)) {
		Py_ssize_t i;
		i = PyNumber_AsSsize_t(item, PyExc_IndexError);

		if (i < 0) {
			i += self->seq_counts;
		}

		if(i >= self->seq_counts){
			PyErr_SetString(PyExc_IndexError, "index out of range");
			return NULL;
		}

		return pyfastx_index_get_seq_by_id(self->index, i+1);
		
	} else if (PyUnicode_CheckExact(item)) {
		char *key = PyUnicode_AsUTF8(item);

		return pyfastx_index_get_seq_by_name(self->index, key);

	} else {
		return NULL;
	}
}

int pyfastx_fasta_length(pyfastx_Fasta *self){
	return self->seq_counts;
}

int pyfastx_fasta_contains(pyfastx_Fasta *self, PyObject *key){
	sqlite3_stmt *stmt;
	
	char *name = PyUnicode_AsUTF8(key);

	sqlite3_prepare_v2(self->index->index_db, "SELECT * FROM seq WHERE chrom=? LIMIT 1;", -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, name, -1, NULL);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		return 0;
	}

	return 1;
}

PyObject *pyfastx_fasta_count(pyfastx_Fasta *self, PyObject *args){
	uint32_t l;
	uint32_t c;
	sqlite3_stmt *stmt;

	if (!PyArg_ParseTuple(args, "I", &l)) {
		return NULL;
	}

	const char *sql = "SELECT COUNT(*) FROM seq WHERE slen>=?";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, l);
	if (sqlite3_step(stmt) == SQLITE_ROW) {
		c = sqlite3_column_int(stmt, 0);
		return Py_BuildValue("I", c);
	}

	Py_RETURN_NONE;
}

PyObject *pyfastx_fasta_nl(pyfastx_Fasta *self, PyObject *args){
	sqlite3_stmt *stmt;
	uint16_t p = 50;
	float half_size;
	uint64_t temp_size = 0;
	uint32_t i = 0;
	uint32_t j = 0;

	if (!PyArg_ParseTuple(args, "|i", &p)) {
		return NULL;
	}

	if (p <= 0 && p >= 100){
		return PyErr_Format(PyExc_ValueError, "the value must between 0 and 100");
	}

	half_size = p/100.0 * self->seq_length;

	const char *sql = "SELECT slen FROM seq ORDER BY slen DESC";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		i++;
		j = sqlite3_column_int(stmt, 0);
		temp_size += j;
		if (temp_size >= half_size) {
			return Py_BuildValue("II", j, i);
		}
	}
	Py_RETURN_NONE;
}

PyObject *pyfastx_fasta_longest(pyfastx_Fasta *self, void* closure){
	sqlite3_stmt *stmt;
	const char *name;
	uint32_t len;
	const char *sql = "SELECT chrom,MAX(slen) FROM seq LIMIT 1";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);

	if (sqlite3_step(stmt) == SQLITE_ROW) {
		name = (const char *)sqlite3_column_text(stmt, 0);
		len = sqlite3_column_int(stmt, 1);
		return Py_BuildValue("sI", name, len);
	}

	Py_RETURN_NONE;
}

PyObject *pyfastx_fasta_shortest(pyfastx_Fasta *self, void* closure){
	sqlite3_stmt *stmt;
	const char *name;
	uint32_t len;
	const char *sql = "SELECT chrom,MIN(slen) FROM seq LIMIT 1";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);

	if (sqlite3_step(stmt) == SQLITE_ROW) {
		name = (const char *)sqlite3_column_text(stmt, 0);
		len = sqlite3_column_int(stmt, 1);
		return Py_BuildValue("sI", name, len);
	}

	Py_RETURN_NONE;
}

PyObject *pyfastx_fasta_mean(pyfastx_Fasta *self, void* closure){
	sqlite3_stmt *stmt;
	double len;
	const char *sql = "SELECT AVG(slen) FROM seq LIMIT 1";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);

	if (sqlite3_step(stmt) == SQLITE_ROW) {
		len = sqlite3_column_double(stmt, 0);
		return Py_BuildValue("d", len);
	}

	Py_RETURN_NONE;
}

PyObject *pyfastx_fasta_median(pyfastx_Fasta *self, void* closure){
	sqlite3_stmt *stmt;
	double m;
	const char *sql;
	if (self->seq_counts % 2 == 0) {
		sql = "SELECT AVG(slen) FROM (SELECT slen FROM seq ORDER BY slen LIMIT ?,2) LIMIT 1";
	} else {
		sql = "SELECT slen FROM seq ORDER BY slen LIMIT ?,1";
	}
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, (self->seq_counts - 1)/2);
	if(sqlite3_step(stmt) == SQLITE_ROW){
		m = sqlite3_column_double(stmt, 0);
		return Py_BuildValue("d", m);
	}

	Py_RETURN_NONE;
}

PyObject *pyfastx_fasta_format(pyfastx_Fasta *self, void* closure) {
	if (self->index->gzip_format) {
		Py_RETURN_TRUE;
	}

	Py_RETURN_FALSE;
}

void pyfastx_fasta_calc_composition(pyfastx_Fasta *self) {
	sqlite3_stmt *stmt;
	const char *sql = "SELECT * FROM comp LIMIT 1";

	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	if (sqlite3_step(stmt) == SQLITE_ROW) {
		return;
	}

	if (sqlite3_exec(self->index->index_db, "BEGIN TRANSACTION;", NULL, NULL, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index->index_db));
		return;
	}

	
	sql = "INSERT INTO comp VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);

	//reading file for kseq
	kstream_t* ks;
	
	//read for line
	kstring_t line = {0, 0, 0};

	//ascii char statistics
	uint32_t seq_comp[128] = {0};
	uint64_t fa_comp[26] = {0};

	uint32_t i;
	uint16_t c;
	uint16_t j;

	Py_BEGIN_ALLOW_THREADS

	ks = ks_init(self->index->gzfd);

	while (ks_getuntil(ks, '\n', &line, 0) >= 0) {
		if (line.s[0] == 62) {
			if (sum_array(seq_comp, 128) > 0) {
				j = 0;
				sqlite3_bind_null(stmt, ++j);
				for (c = 65; c <= 90; c++) {
					sqlite3_bind_int(stmt, ++j, seq_comp[c]+seq_comp[c+32]);
					fa_comp[c-65] += seq_comp[c]+seq_comp[c+32];
				}
				sqlite3_step(stmt);
				sqlite3_reset(stmt);
			}
			memset(seq_comp, 0, sizeof(seq_comp));
			continue;
		}

		for (i = 0; i < line.l; i++) {
			c = line.s[i];
			if (c == 10 || c == 13) {
				continue;
			}
			++seq_comp[c];
		}
	}

	if (sum_array(seq_comp, 128) > 0) {
		j = 0;
		sqlite3_bind_null(stmt, ++j);
		for (c = 65; c <= 90; c++) {
			sqlite3_bind_int(stmt, ++j, seq_comp[c]+seq_comp[c+32]);
			fa_comp[c-65] += seq_comp[c]+seq_comp[c+32];
		}
		sqlite3_step(stmt);
		sqlite3_reset(stmt);
	}

	//write total bases to db
	sqlite3_bind_null(stmt, 1);
	for (j = 0; j < 26; j++) {
		sqlite3_bind_int64(stmt, j+2, fa_comp[j]);
	}
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
	sqlite3_exec(self->index->index_db, "COMMIT;", NULL, NULL, NULL);

	ks_destroy(ks);
	free(line.s);

	Py_END_ALLOW_THREADS
}

PyObject *pyfastx_fasta_gc_content(pyfastx_Fasta *self, void* closure) {
	sqlite3_stmt *stmt;
	int64_t a, c, g, t;
	pyfastx_fasta_calc_composition(self);
	const char *sql = "SELECT a, c, g, t FROM comp ORDER BY ID DESC LIMIT 1";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	
	if (sqlite3_step(stmt) != SQLITE_ROW) {
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index->index_db));
		return NULL;
	}

	a = sqlite3_column_int64(stmt, 0);
	c = sqlite3_column_int64(stmt, 1);
	g = sqlite3_column_int64(stmt, 2);
	t = sqlite3_column_int64(stmt, 3);

	return Py_BuildValue("f", (float)(g+c)/(a+c+g+t)*100);
}

PyObject *pyfastx_fasta_gc_skew(pyfastx_Fasta *self, void* closure) {
	sqlite3_stmt *stmt;
	int64_t c, g;
	pyfastx_fasta_calc_composition(self);
	const char *sql = "SELECT c, g FROM comp ORDER BY ID DESC LIMIT 1";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	
	if (sqlite3_step(stmt) != SQLITE_ROW) {
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index->index_db));
		return NULL;
	}
	
	c = sqlite3_column_int64(stmt, 0);
	g = sqlite3_column_int64(stmt, 1);
	
	return Py_BuildValue("f", (float)(g-c)/(g+c));
}

PyObject *pyfastx_fasta_composition(pyfastx_Fasta *self, void* closure) {
	sqlite3_stmt *stmt;
	int i;
	int64_t c;
	pyfastx_fasta_calc_composition(self);
	const char *sql = "SELECT * FROM comp ORDER BY ID DESC LIMIT 1";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, self->seq_counts+1);
	
	if (sqlite3_step(stmt) != SQLITE_ROW) {
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index->index_db));
		return NULL;
	}

	PyObject *d = PyDict_New();

	for (i = 1; i < 27; i++) {
		c = sqlite3_column_int64(stmt, i);
		if (c > 0) {
			PyDict_SetItem(d, Py_BuildValue("s", int_to_str(i+64)), Py_BuildValue("l", c));
		}
	}

	return d;
}

//support for guess sequence type according to IUPAC codes
//https://www.bioinformatics.org/sms/iupac.html
PyObject *pyfastx_fasta_guess_type(pyfastx_Fasta *self, void* closure) {
	sqlite3_stmt *stmt;
	int i, j;
	int64_t c;
	char *alphabets;
	char *ret;

	pyfastx_fasta_calc_composition(self);

	const char *sql = "SELECT * FROM comp ORDER BY ID DESC LIMIT 1";
	sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);

	if (sqlite3_step(stmt) != SQLITE_ROW) {
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index->index_db));
		return NULL;
	}

	alphabets = (char *)malloc(26);
	j = 0;
	for (i = 1; i < 27; i++) {
		c = sqlite3_column_int64(stmt, i);
		if (c > 0) {
			alphabets[j++] = i+64;
		}
	}
	alphabets[j] = '\0';

	if (is_subset("ACGTN", alphabets) || is_subset("ABCDGHKMNRSTVWY", alphabets)) {
		ret = "DNA";
	} else if (is_subset("ACGUN", alphabets) || is_subset("ABCDGHKMNRSUVWY", alphabets)) {
		ret = "RNA";
	} else if (is_subset("ACDEFGHIKLMNPQRSTVWY", alphabets)) {
		ret = "protein";
	} else {
		ret = "unknown";
	}

	return Py_BuildValue("s", ret);
}

static PyGetSetDef pyfastx_fasta_getsets[] = {
	{"longest", (getter)pyfastx_fasta_longest, NULL, NULL, NULL},
	{"shortest", (getter)pyfastx_fasta_shortest, NULL, NULL, NULL},
	{"mean", (getter)pyfastx_fasta_mean, NULL, NULL, NULL},
	{"median", (getter)pyfastx_fasta_median, NULL, NULL, NULL},
	{"is_gzip", (getter)pyfastx_fasta_format, NULL, NULL, NULL},
	{"gc_content", (getter)pyfastx_fasta_gc_content, NULL, NULL, NULL},
	{"gc_skew", (getter)pyfastx_fasta_gc_skew, NULL, NULL, NULL},
	{"composition", (getter)pyfastx_fasta_composition, NULL, NULL, NULL},
	{"type", (getter)pyfastx_fasta_guess_type, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef pyfastx_fasta_members[] = {
	{"file_name", T_STRING, offsetof(pyfastx_Fasta, file_name), READONLY},
	{"size", T_LONG, offsetof(pyfastx_Fasta, seq_length), READONLY},
	//{"count", T_INT, offsetof(pyfastx_Fasta, seq_counts), READONLY},
	//{"gc_content", T_FLOAT, offsetof(pyfastx_Fasta, gc_content), READONLY},
	//{"gc_skew", T_FLOAT, offsetof(pyfastx_Fasta, gc_skew), READONLY},
	//{"composition", T_OBJECT, offsetof(pyfastx_Fasta, composition), READONLY},
	{NULL}
};

static PyMethodDef pyfastx_fasta_methods[] = {
	{"build_index", (PyCFunction)pyfastx_fasta_build_index, METH_VARARGS},
	{"rebuild_index", (PyCFunction)pyfastx_fasta_rebuild_index, METH_VARARGS},
	{"fetch", (PyCFunction)pyfastx_fasta_fetch, METH_VARARGS|METH_KEYWORDS},
	{"count", (PyCFunction)pyfastx_fasta_count, METH_VARARGS},
	{"keys", (PyCFunction)pyfastx_fasta_keys, METH_VARARGS},
	{"nl", (PyCFunction)pyfastx_fasta_nl, METH_VARARGS},
	//{"test", (PyCFunction)test, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

//as a list
static PySequenceMethods seq_methods = {
	0, /*sq_length*/
	0, /*sq_concat*/
	0, /*sq_repeat*/
	0, /*sq_item*/
	0, /*sq_slice */
	0, /*sq_ass_item*/
	0, /*sq_ass_splice*/
	(objobjproc)pyfastx_fasta_contains, /*sq_contains*/
	0, /*sq_inplace_concat*/
	0, /*sq_inplace_repeat*/
};

static PyMappingMethods pyfastx_fasta_as_mapping = {
	(lenfunc)pyfastx_fasta_length,
	(binaryfunc)pyfastx_fasta_subscript,
	0,
};

PyTypeObject pyfastx_FastaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Fasta",                        /* tp_name */
    sizeof(pyfastx_Fasta),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)pyfastx_fasta_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)pyfastx_fasta_repr,                              /* tp_repr */
    0,                              /* tp_as_number */
    &seq_methods,                   /* tp_as_sequence */
    &pyfastx_fasta_as_mapping,                   /* tp_as_mapping */
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
    (getiterfunc)pyfastx_fasta_iter,     /* tp_iter */
    (iternextfunc)pyfastx_fasta_next,    /* tp_iternext */
    pyfastx_fasta_methods,          /* tp_methods */
    pyfastx_fasta_members,          /* tp_members */
    pyfastx_fasta_getsets,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    pyfastx_fasta_new,              /* tp_new */
};
