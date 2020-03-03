#include "fasta.h"
#include "util.h"
#include "identifier.h"
#include "structmember.h"

/*calculate fasta attributes including sequence count, length,
composition (ATGCN count) and GC content
*/
void pyfastx_calc_fasta_attrs(pyfastx_Fasta *self){
	sqlite3_stmt *stmt;
	int ret;
	
	//sequence count
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, "SELECT * FROM stat LIMIT 1;", -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			self->seq_counts = sqlite3_column_int(stmt, 0);
			self->seq_length = sqlite3_column_int64(stmt, 1);
		);
	} else {
		PyErr_SetString(PyExc_RuntimeError, "get seq count and length error");
		return;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
}

PyObject *pyfastx_fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	//fasta file path
	char *file_name;
	int file_len;

	//bool value for uppercase sequence
	int uppercase = 1;

	//build index or not
	int build_index = 1;

	//just keep index into memory, do not generate index file
	int memory_index = 0;

	//calculate the composition of sequence
	int full_index = 0;

	//key function for seperating name
	PyObject *key_func = Py_None;

	pyfastx_Fasta *obj;

	//paramters for fasta object construction
	static char* keywords[] = {"file_name", "uppercase", "build_index", "full_index", "memory_index", "key_func", NULL};
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s#|iiiiO", keywords, &file_name, &file_len, &uppercase, &build_index, &full_index, &memory_index, &key_func)){
		return NULL;
	}

	if ((key_func != Py_None) && !PyCallable_Check(key_func)) {
		PyErr_SetString(PyExc_TypeError, "key_func must be a callable function");
		return NULL;
	}

	//check input sequence file is whether exists
	if (!file_exists(file_name)) {
		PyErr_Format(PyExc_FileExistsError, "input fasta file %s does not exists", file_name);
		return NULL;
	}

	//create Fasta class
	obj = (pyfastx_Fasta *)type->tp_alloc(type, 0);
	if (!obj){
		return NULL;
	}
	
	//initial sequence file name
	obj->file_name = (char *)malloc(file_len + 1);
	strcpy(obj->file_name, file_name);

	obj->uppercase = uppercase;
	obj->has_index = build_index;

	//create index
	obj->index = pyfastx_init_index(obj->file_name, file_len, uppercase, memory_index, key_func);

	//initial iterator stmt
	obj->iter_stmt = NULL;

	//check is correct fasta format
	if (!fasta_validator(obj->index->gzfd)) {
		PyErr_Format(PyExc_RuntimeError, "%s is not plain or gzip compressed fasta formatted file", file_name);
		return NULL;
	}

	//if build_index is True
	if (build_index) {
		pyfastx_build_index(obj->index);
		pyfastx_calc_fasta_attrs(obj);

		if (full_index) {
			pyfastx_fasta_calc_composition(obj);
		}
	}
	
	return (PyObject *)obj;
}

void pyfastx_fasta_dealloc(pyfastx_Fasta *self){
	pyfastx_index_free(self->index);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *pyfastx_fasta_iter(pyfastx_Fasta *self){
	pyfastx_rewind_index(self->index);

	if (self->has_index) {
		//self->iter_id = 0;
		PYFASTX_SQLITE_CALL(
			sqlite3_finalize(self->iter_stmt);
			self->iter_stmt = NULL;
			sqlite3_prepare_v2(self->index->index_db, "SELECT * FROM seq", -1, &self->iter_stmt, NULL);
		);
	}

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fasta_repr(pyfastx_Fasta *self){
	return PyUnicode_FromFormat("<Fasta> %s contains %d sequences", self->file_name, self->seq_counts);
}

PyObject *pyfastx_fasta_next(pyfastx_Fasta *self){
	if (self->has_index) {
		//++self->iter_id;

		//if (self->iter_id <= self->seq_counts) {
		//	return pyfastx_index_get_seq_by_id(self->index, self->iter_id);
		//}

		int ret;
		PYFASTX_SQLITE_CALL(ret = sqlite3_step(self->iter_stmt));
		if (ret == SQLITE_ROW) {
			return pyfastx_index_make_seq(self->index, self->iter_stmt);
		}

	} else {
		return pyfastx_get_next_seq(self->index);
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(self->iter_stmt));

	return NULL;
}

PyObject *pyfastx_fasta_build_index(pyfastx_Fasta *self){
	if (!self->index->index_db) {
		pyfastx_build_index(self->index);
		pyfastx_calc_fasta_attrs(self);
	}

	Py_RETURN_TRUE;
}

/*PyObject *pyfastx_fasta_rebuild_index(pyfastx_Fasta *self){
	if (self->index->index_db) {
		PYFASTX_SQLITE_CALL(sqlite3_close(self->index->index_db));
		self->index->index_db = 0;
	}
	
	if (file_exists(self->index->index_file)) {
		remove(self->index->index_file);
	}

	pyfastx_build_index(self->index);
	pyfastx_calc_fasta_attrs(self);
	Py_RETURN_TRUE;
}*/

PyObject *pyfastx_fasta_fetch(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	static char* keywords[] = {"name", "intervals", "strand", NULL};

	char *name;
	char *seq;
	PyObject *intervals;
	uint64_t start;
	uint64_t end;
	int strand = '+';
	PyObject *item;
	Py_ssize_t size;
	sqlite3_stmt *stmt;
	const char* sql;
	uint32_t chrom;
	uint32_t seq_len;
	char* sub_seq;
	int ret;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "sO|C", keywords, &name, &intervals, &strand)){
		return NULL;
	}

	if(!PyTuple_Check(intervals) && !PyList_Check(intervals)){
		PyErr_SetString(PyExc_ValueError, "intervals must be list or tuple");
		return NULL;
	}

	if (PyList_Check(intervals)) {
		intervals = PyList_AsTuple(intervals);
	}

	item = PyTuple_GetItem(intervals, 0);
	size = PyTuple_Size(intervals);
	
	//select sql statement, chrom indicates seq name or chromomsome
	sql = "SELECT ID FROM seq WHERE chrom=? LIMIT 1;";

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
	);
	
	if (ret == SQLITE_ROW){
		PYFASTX_SQLITE_CALL(chrom = sqlite3_column_int(stmt, 0));
	} else {
		PyErr_Format(PyExc_NameError, "Sequence %s does not exists", name);
		return NULL;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	seq = pyfastx_index_get_full_seq(self->index, chrom);

	if (PyLong_Check(item)) {
		if (size != 2) {
			PyErr_SetString(PyExc_ValueError, "list or tuple should include only start and end");
			return NULL;
		}

		start = PyLong_AsLong(item);
		item = PyTuple_GetItem(intervals, 1);
		end = PyLong_AsLong(item);

		if (start > end) {
			PyErr_SetString(PyExc_ValueError, "start position should less than end position");
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

		for (i=0; i<size; i++) {
			item = PyTuple_GetItem(intervals, i);
			
			if (PyList_Check(item)) {
				item = PyList_AsTuple(item);
			}

			start = PyLong_AsLong(PyTuple_GetItem(item, 0));
			end = PyLong_AsLong(PyTuple_GetItem(item, 1));
			seq_len = end - start + 1;

			if (start > end) {
				PyErr_SetString(PyExc_ValueError, "start position should less than end position");
				return NULL;
			}

			memcpy(sub_seq+j, seq+start-1, seq_len);
			j += seq_len;
		}
		sub_seq[j] = '\0';
	}

	if (strand == '-') {
		reverse_complement_seq(sub_seq);
	}
	
	return Py_BuildValue("s", sub_seq);
}

PyObject *pyfastx_fasta_keys(pyfastx_Fasta *self) {
	//pyfastx_Identifier *ids = PyObject_New(pyfastx_Identifier, &pyfastx_IdentifierType);
	pyfastx_Identifier *ids = (pyfastx_Identifier *)PyObject_CallObject((PyObject *)&pyfastx_IdentifierType, NULL);
	
	if (!ids) {
		return NULL;
	}

	ids->index_db = self->index->index_db;
	ids->stmt = NULL;
	ids->seq_counts = self->seq_counts;
	ids->sort = 0;
	ids->order = 0;
	ids->filter = NULL;
	ids->temp_filter = NULL;

	//Py_INCREF(ids);
	return (PyObject *)ids;
}

PyObject *pyfastx_fasta_subscript(pyfastx_Fasta *self, PyObject *item){
	
	if (PyIndex_Check(item)) {
		Py_ssize_t i;
		i = PyNumber_AsSsize_t(item, PyExc_IndexError);

		if (i < 0) {
			i += self->seq_counts;
		}

		if (i >= self->seq_counts) {
			PyErr_SetString(PyExc_IndexError, "index out of range");
			return NULL;
		}

		return pyfastx_index_get_seq_by_id(self->index, i+1);
		
	} else if (PyUnicode_CheckExact(item)) {
		return pyfastx_index_get_seq_by_name(self->index, PyUnicode_AsUTF8(item));

	} else {
		PyErr_SetString(PyExc_KeyError, "the key must be index number or sequence name");
		return NULL;
	}
}

int pyfastx_fasta_length(pyfastx_Fasta *self){
	return self->seq_counts;
}

int pyfastx_fasta_contains(pyfastx_Fasta *self, PyObject *key){
	sqlite3_stmt *stmt;
	int ret;
	char *name;

	if (!PyUnicode_CheckExact(key)) {
		return 0;
	}
	
	name = PyUnicode_AsUTF8(key);

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, "SELECT * FROM seq WHERE chrom=? LIMIT 1;", -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	return ret == SQLITE_ROW ? 1 : 0;
}

PyObject *pyfastx_fasta_count(pyfastx_Fasta *self, PyObject *args){
	int l;
	sqlite3_stmt *stmt;
	int ret;
	int c;

	const char *sql;
	if (!PyArg_ParseTuple(args, "i", &l)) {
		return NULL;
	}

	sql = "SELECT COUNT(*) FROM seq WHERE slen>=?;";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int(stmt, 1, l);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(c = sqlite3_column_int(stmt, 0));
	} else {
		c = 0;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return Py_BuildValue("i", c);
}

PyObject *pyfastx_fasta_nl(pyfastx_Fasta *self, PyObject *args){
	sqlite3_stmt *stmt;
	int p = 50;
	float half_size;
	uint64_t temp_size = 0;
	uint32_t i = 0;
	uint32_t j = 0;
	const char *sql;
	int ret;

	if (!PyArg_ParseTuple(args, "|i", &p)) {
		return NULL;
	}

	if (p < 0 || p > 100){
		PyErr_SetString(PyExc_ValueError, "the value must between 0 and 100");
		return NULL;
	}

	half_size = p/100.0 * self->seq_length;

	sql = "SELECT slen FROM seq ORDER BY slen DESC";
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL));
	for(;;){
		PYFASTX_SQLITE_CALL(ret=sqlite3_step(stmt));
		if (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(j = sqlite3_column_int(stmt, 0));
			i++;
			temp_size += j;
			if (temp_size >= half_size)
				break;
		} else {
			break;
		}
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return Py_BuildValue("II", j, i);
}

PyObject *pyfastx_fasta_longest(pyfastx_Fasta *self, void* closure){
	sqlite3_stmt *stmt;
	int ret;
	uint32_t chrom;
	const char *sql;

	sql = "SELECT ID,MAX(slen) FROM seq LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			chrom = sqlite3_column_int(stmt, 0);
			sqlite3_finalize(stmt);
		);

		return pyfastx_index_get_seq_by_id(self->index, chrom);
	}
	
	PyErr_SetString(PyExc_RuntimeError, "not found longest sequence");
	return NULL;
}

PyObject *pyfastx_fasta_shortest(pyfastx_Fasta *self, void* closure){
	sqlite3_stmt *stmt;
	int ret;
	uint32_t chrom;
	const char *sql;
	
	sql = "SELECT ID,MIN(slen) FROM seq LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			chrom = sqlite3_column_int(stmt, 0);
			sqlite3_finalize(stmt);
		);

		return pyfastx_index_get_seq_by_id(self->index, chrom);
	}

	PyErr_SetString(PyExc_RuntimeError, "not found shortest sequence");
	return NULL;
}

PyObject *pyfastx_fasta_mean(pyfastx_Fasta *self, void* closure){
	sqlite3_stmt *stmt;
	int ret;
	double len;
	const char *sql;

	sql = "SELECT AVG(slen) FROM seq LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			len = sqlite3_column_double(stmt, 0);
			sqlite3_finalize(stmt);
		);
		return Py_BuildValue("d", len);
	}

	PyErr_SetString(PyExc_RuntimeError, "can not calculate average length");
	return NULL;
}

PyObject *pyfastx_fasta_median(pyfastx_Fasta *self, void* closure){
	sqlite3_stmt *stmt;
	int ret;
	double m;
	const char *sql;

	if (self->seq_counts % 2 == 0) {
		sql = "SELECT AVG(slen) FROM (SELECT slen FROM seq ORDER BY slen LIMIT ?,2) LIMIT 1";
	} else {
		sql = "SELECT slen FROM seq ORDER BY slen LIMIT ?,1";
	}

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int(stmt, 1, (self->seq_counts - 1)/2);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW){
		PYFASTX_SQLITE_CALL(
			m = sqlite3_column_double(stmt, 0);
			sqlite3_finalize(stmt);
		);
		return Py_BuildValue("d", m);
	}

	PyErr_SetString(PyExc_RuntimeError, "can not calculate median length");
	return NULL;
}

PyObject *pyfastx_fasta_format(pyfastx_Fasta *self, void* closure) {
	if (self->index->gzip_format) {
		Py_RETURN_TRUE;
	}

	Py_RETURN_FALSE;
}

void pyfastx_fasta_calc_composition(pyfastx_Fasta *self) {
	sqlite3_stmt *stmt;
	int ret;
	const char *sql;

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
	
	sql = "SELECT * FROM comp LIMIT 1";

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);
	
	if (ret == SQLITE_ROW)
		return;

	stmt = NULL;

	PYFASTX_SQLITE_CALL(sqlite3_exec(self->index->index_db, "BEGIN TRANSACTION;", NULL, NULL, NULL));

	sql = "INSERT INTO comp VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL));
	
	gzrewind(self->index->gzfd);
	ks = ks_init(self->index->gzfd);
	
	Py_BEGIN_ALLOW_THREADS
	
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

	Py_END_ALLOW_THREADS
	
	ks_destroy(ks);
	free(line.s);
}

PyObject *pyfastx_fasta_gc_content(pyfastx_Fasta *self, void* closure) {
	sqlite3_stmt *stmt;
	int ret;
	int64_t a, c, g, t;
	const char *sql;

	pyfastx_fasta_calc_composition(self);
	sql = "SELECT a, c, g, t FROM comp ORDER BY ID DESC LIMIT 1";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			a = sqlite3_column_int64(stmt, 0);
			c = sqlite3_column_int64(stmt, 1);
			g = sqlite3_column_int64(stmt, 2);
			t = sqlite3_column_int64(stmt, 3);
			sqlite3_finalize(stmt);
		);
		return Py_BuildValue("f", (float)(g+c)/(a+c+g+t)*100);
	
	} else {
		PyErr_SetString(PyExc_RuntimeError, "can not calculate gc content");
		return NULL;
	}
}

PyObject *pyfastx_fasta_gc_skew(pyfastx_Fasta *self, void* closure) {
	sqlite3_stmt *stmt;
	int ret;
	int64_t c, g;
	const char *sql;

	pyfastx_fasta_calc_composition(self);
	sql = "SELECT c, g FROM comp ORDER BY ID DESC LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);
	
	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			c = sqlite3_column_int64(stmt, 0);
			g = sqlite3_column_int64(stmt, 1);
			sqlite3_finalize(stmt);
		);
		return Py_BuildValue("f", (float)(g-c)/(g+c));
	} else {
		PyErr_SetString(PyExc_RuntimeError, "can not calculate gc skew");
		return NULL;
	}
}

PyObject *pyfastx_fasta_composition(pyfastx_Fasta *self, void* closure) {
	sqlite3_stmt *stmt;
	int ret;
	int i;
	int64_t c;
	const char *sql;
	PyObject *d;

	pyfastx_fasta_calc_composition(self);

	//the last row store the sum of the each base
	sql = "SELECT * FROM comp ORDER BY ID DESC LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);
	
	if (ret == SQLITE_ROW) {
		d = PyDict_New();
		for (i = 1; i < 27; i++) {
			PYFASTX_SQLITE_CALL(c = sqlite3_column_int64(stmt, i));
			if (c > 0) {
				PyDict_SetItem(d, Py_BuildValue("C", i+64), Py_BuildValue("l", c));
			}
		}
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		return d;
		
	} else {
		PyErr_SetString(PyExc_RuntimeError, "can not get composition");
		return NULL;
	}
}

//support for guess sequence type according to IUPAC codes
//https://www.bioinformatics.org/sms/iupac.html
PyObject *pyfastx_fasta_guess_type(pyfastx_Fasta *self, void* closure) {
	sqlite3_stmt *stmt;
	int ret;
	int i, j;
	char *alphabets;
	char *retval;
	const char *sql;
	int64_t c;

	pyfastx_fasta_calc_composition(self);

	sql = "SELECT * FROM comp ORDER BY ID DESC LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret != SQLITE_ROW) {
		PyErr_SetString(PyExc_RuntimeError, "can not get sequence type");
		return NULL;
	}

	alphabets = (char *)malloc(26);
	j = 0;
	for (i = 1; i < 27; i++) {
		PYFASTX_SQLITE_CALL(c = sqlite3_column_int64(stmt, i));
		if (c > 0) {
			alphabets[j++] = i+64;
		}
	}
	alphabets[j] = '\0';
	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	if (is_subset("ACGTN", alphabets) || is_subset("ABCDGHKMNRSTVWY-", alphabets)) {
		retval = "DNA";
	} else if (is_subset("ACGUN", alphabets) || is_subset("ABCDGHKMNRSUVWY-", alphabets)) {
		retval = "RNA";
	} else if (is_subset("ACDEFGHIKLMNPQRSTVWY*-", alphabets)) {
		retval = "protein";
	} else {
		retval = "unknown";
	}

	return Py_BuildValue("s", retval);
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
	{"build_index", (PyCFunction)pyfastx_fasta_build_index, METH_NOARGS, NULL},
	//{"rebuild_index", (PyCFunction)pyfastx_fasta_rebuild_index, METH_NOARGS, NULL},
	{"fetch", (PyCFunction)pyfastx_fasta_fetch, METH_VARARGS|METH_KEYWORDS, NULL},
	{"count", (PyCFunction)pyfastx_fasta_count, METH_VARARGS, NULL},
	{"keys", (PyCFunction)pyfastx_fasta_keys, METH_NOARGS, NULL},
	{"nl", (PyCFunction)pyfastx_fasta_nl, METH_VARARGS, NULL},
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
