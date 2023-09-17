#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "fasta.h"
#include "util.h"
#include "fakeys.h"
#include "structmember.h"
#include "sequence.h"
#include "stdint.h"
#include "zlib.h"
#include "kseq.h"
#include "zran.h"
#include "sqlite3.h"

/*calculate fasta attributes including sequence count, length,
composition (ATGCN count) and GC content
*/
void pyfastx_calc_fasta_attrs(pyfastx_Fasta *self){
	int ret;
	sqlite3_stmt *stmt;

	//sequence count
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, "SELECT * FROM stat LIMIT 1", -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			self->seq_counts = sqlite3_column_int64(stmt, 0);
			self->seq_length = sqlite3_column_int64(stmt, 1);
		);
	} else {
		PyErr_SetString(PyExc_RuntimeError, "get seq count and length error");
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
}

PyObject *pyfastx_fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
	//bool value for uppercase sequence
	int uppercase = 0;

	//build index or not
	int build_index = 1;

	//just keep index into memory, do not generate index file
	int memory_index = 0;

	//calculate the composition of sequence
	int full_index = 0;

	//use full name instead of identifier before first whitespace
	int full_name = 0;

	//fasta file path
	PyObject *file_obj;
	PyObject *index_obj = NULL;

	//key function for seperating name
	PyObject *key_func = NULL;

	pyfastx_Fasta *obj;

	//paramters for fasta object construction
	static char* keywords[] = {"file_name", "index_file", "uppercase", "build_index", "full_index", "full_name", "memory_index", "key_func", NULL};
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O|OiiiiiO", keywords, &file_obj, &index_obj, &uppercase, &build_index, &full_index, &full_name, &memory_index, &key_func)){
		return NULL;
	}

	if ((key_func) && !PyCallable_Check(key_func)) {
		PyErr_SetString(PyExc_TypeError, "key_func must be a callable function");
		return NULL;
	}

	//check input sequence file is whether exists
	//file_name = (char *)PyUnicode_AsUTF8AndSize(file_obj, &file_len);

	/*if (!file_name) {
		PyErr_Format(PyExc_ValueError, "the input file name is not a right string");
		return NULL;
	}*/

	if (!file_exists(file_obj)) {
		PyErr_Format(PyExc_FileExistsError, "the input fasta file %U does not exists", file_obj);
		return NULL;
	}

	//create Fasta class
	obj = (pyfastx_Fasta *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	//initial sequence file name
	obj->file_obj = Py_NewRef(file_obj);

	obj->uppercase = uppercase;
	obj->has_index = build_index;

	//create index

	obj->index = pyfastx_init_index((PyObject *)obj, file_obj, index_obj, uppercase, full_name, memory_index, key_func);
	
	//iter function
	obj->func = pyfastx_index_next_null;

	//check is correct fasta format
	if (!fasta_validator(obj->index->gzfd)) {
		PyErr_Format(PyExc_RuntimeError, "%U is not plain or gzip compressed fasta formatted file", file_obj);
		return NULL;
	}

	//if build_index is True
	if (build_index) {
		pyfastx_build_index(obj->index);
		pyfastx_calc_fasta_attrs(obj);

		if (full_index) {
			pyfastx_fasta_calc_composition(obj);
		}

		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(obj->index->index_db, "SELECT * FROM seq WHERE chrom=? LIMIT 1;", -1, &obj->index->seq_stmt, NULL);
			sqlite3_prepare_v2(obj->index->index_db, "SELECT * FROM seq WHERE ID=? LIMIT 1;", -1, &obj->index->uid_stmt, NULL);
		);
	}

	return (PyObject *)obj;
}

void pyfastx_fasta_dealloc(pyfastx_Fasta *self){
	//free(self->file_name);
	pyfastx_index_free(self->index);
	Py_DECREF(self->file_obj);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *pyfastx_fasta_repr(pyfastx_Fasta *self) {
	if (self->has_index) {
		return PyUnicode_FromFormat("<Fasta> %U contains %zd sequences", self->file_obj, self->seq_counts);
	} else {
		return PyUnicode_FromFormat("<Fasta> %U", self->file_obj);
	}
}

PyObject *pyfastx_fasta_iter(pyfastx_Fasta *self){
	pyfastx_rewind_index(self->index);

	if (self->has_index) {
		self->index->iterating = 1;
		PYFASTX_SQLITE_CALL(
			sqlite3_finalize(self->index->iter_stmt);
			self->index->iter_stmt = NULL;
			sqlite3_prepare_v2(self->index->index_db, "SELECT * FROM seq", -1, &self->index->iter_stmt, NULL);
		);

		self->func = pyfastx_index_next_with_index_seq;
	} else {
		if (self->index->uppercase && self->index->full_name) {
			self->func = pyfastx_index_next_full_name_upper_seq;
		} else if (self->index->uppercase) {
			self->func = pyfastx_index_next_upper_seq;
		} else if (self->index->full_name) {
			self->func = pyfastx_index_next_full_name_seq;
		} else {
			self->func = pyfastx_index_next_seq;
		}
	}

	Py_INCREF(self);
	return (PyObject *)self;
}

PyObject *pyfastx_fasta_next(pyfastx_Fasta *self) {
	return self->func(self->index);
}

PyObject *pyfastx_fasta_build_index(pyfastx_Fasta *self){
	if (!self->index->index_db) {
		pyfastx_build_index(self->index);
		pyfastx_calc_fasta_attrs(self);
		self->has_index = 1;
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

PyObject * pyfastx_fasta_slice_from_cache(pyfastx_Fasta *self, Py_ssize_t start, Py_ssize_t end, int flank) {
	char *left;
	char *right;

	Py_ssize_t slice_start;
	Py_ssize_t slice_len;
	
	PyObject *ret;

	slice_start = start - flank - 1;
	if (slice_start < 0) {
		slice_len = flank + slice_start;
		slice_start = 0;
	} else {
		slice_len = flank;
	}

	if (slice_len > 0) {
		left = (char *)malloc(slice_len + 1);
		memcpy(left, self->index->cache_seq.s+slice_start, slice_len);
		left[slice_len] = '\0';
	} else {
		left = (char *)malloc(1);
		left[0] = '\0';
	}

	slice_start = end;
	if ((end + flank) > self->index->cache_seq.l) {
		slice_len = self->index->cache_seq.l - end;
	} else {
		slice_len = flank;
	}

	if (slice_len > 0) {
		right = (char *)malloc(slice_len + 1);
		memcpy(right, self->index->cache_seq.s+slice_start, slice_len);
		right[slice_len] = '\0';
	} else {
		right = (char *)malloc(1);
		right[0] = '\0';
	}

	ret = Py_BuildValue("ss", left, right);
	free(left);
	free(right);

	return ret;
}

void pyfastx_fasta_seq_info(pyfastx_Fasta *self, char *name, Py_ssize_t *chrom, Py_ssize_t *offset, Py_ssize_t *bytes, Py_ssize_t *slen, Py_ssize_t *llen, int *elen, int *normal) {
	int ret;
	sqlite3_stmt *stmt;

	const char *sql = "SELECT ID,boff,blen,slen,llen,elen,norm FROM seq WHERE chrom=? LIMIT 1;";

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
	);
	
	if (ret == SQLITE_ROW){
		PYFASTX_SQLITE_CALL(
			*chrom = sqlite3_column_int64(stmt, 0);
			*offset = sqlite3_column_int64(stmt, 1);
			*bytes = sqlite3_column_int64(stmt, 2);
			*slen = sqlite3_column_int64(stmt, 3);
			*llen = sqlite3_column_int64(stmt, 4);
			*elen = sqlite3_column_int(stmt, 5);
			*normal = sqlite3_column_int(stmt, 6);
		);
	} else {
		PyErr_Format(PyExc_NameError, "sequence %s does not exists", name);
	}
	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
}

void pyfastx_fasta_cache_full(pyfastx_Fasta *self, Py_ssize_t chrom, Py_ssize_t offset, Py_ssize_t bytes) {
	pyfastx_index_fill_cache(self->index, offset, bytes);
	self->index->cache_chrom = chrom;
	self->index->cache_start = 1;
	self->index->cache_end = self->index->cache_seq.l;
	self->index->cache_full = 1;
}

char *pyfastx_fasta_slice_seq(pyfastx_Fasta *self, Py_ssize_t offset, Py_ssize_t bytelen, Py_ssize_t line_len, int end_len, Py_ssize_t slice_start, Py_ssize_t slice_stop) {
	char *ret;

	Py_ssize_t before_sline;
	Py_ssize_t before_eline;
	Py_ssize_t cross_line;

	if (slice_stop >= slice_start) {
		before_sline = slice_start/(line_len - end_len);
		before_eline = slice_stop/(line_len - end_len);
		cross_line = before_eline - before_sline;
		offset = offset + slice_start + end_len*before_sline;
		bytelen = slice_stop - slice_start + cross_line*end_len;

		ret = (char *)malloc(bytelen + 1);
		pyfastx_index_random_read(self->index, ret, offset, bytelen);
		if (self->index->uppercase) {
			remove_space_uppercase(ret, bytelen);
		} else {
			remove_space(ret, bytelen);
		}
	} else {
		ret = (char *)malloc(1);
		ret[0] = '\0';
	}

	return ret;
}

PyObject *pyfastx_fasta_flank(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs) {
	int end_len;
	int normal;
	int use_cache = 0;
	int flank_len = 50;

	char *name;
	char *left;
	char *right;

	Py_ssize_t start;
	Py_ssize_t end;
	Py_ssize_t slice_start;
	Py_ssize_t slice_stop;
	Py_ssize_t chrom;
	Py_ssize_t offset;
	Py_ssize_t bytes;
	Py_ssize_t seq_len;
	Py_ssize_t line_len;

	PyObject *ret;

	static char *keywords[] = {"chrom", "start", "end", "flank_length", "use_cache", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "snn|ii", keywords, &name, &start, &end, &flank_len, &use_cache)) {
		return NULL;
	}

	if (self->index->cache_name.s && strcmp(self->index->cache_name.s, name) == 0 && self->index->cache_full) {
		return pyfastx_fasta_slice_from_cache(self, start, end, flank_len);
	}

	pyfastx_fasta_seq_info(self, name, &chrom, &offset, &bytes, &seq_len, &line_len, &end_len, &normal);

	if (use_cache || !normal) {
		pyfastx_fasta_cache_full(self, chrom, offset, bytes);
		ret = pyfastx_fasta_slice_from_cache(self, start, end, flank_len);
	} else {
		slice_start = start - flank_len - 1;
		if (slice_start < 0) {
			slice_start = 0;
		}
		slice_stop = start - 1;

		//slice_start 0-base, slice_stop 1-base
		left = pyfastx_fasta_slice_seq(self, offset, bytes, line_len, end_len, slice_start, slice_stop);

		slice_stop = end + flank_len;
		if (slice_stop > seq_len) {
			slice_stop = seq_len;
		}
		slice_start = end;

		right = pyfastx_fasta_slice_seq(self, offset, bytes, line_len, end_len, slice_start, slice_stop);
		ret = Py_BuildValue("ss", left, right);
		free(left);
		free(right);
	}

	return ret;
}

PyObject *pyfastx_fasta_fetch(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
	static char* keywords[] = {"chrom", "intervals", "strand", NULL};

	char *name;
	char *seq;
	char* sub_seq;
	const char* sql;

	int ret;
	int strand = '+';

	sqlite3_stmt *stmt;
	
	Py_ssize_t start;
	Py_ssize_t end;
	Py_ssize_t size;
	Py_ssize_t chrom;
	Py_ssize_t seq_len;
	Py_ssize_t offset;
	Py_ssize_t bytes;
	PyObject *intervals;
	PyObject *item;

	
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
	if (self->index->cache_name.s && strcmp(self->index->cache_name.s, name) == 0 && self->index->cache_full) {
		seq = self->index->cache_seq.s;
	} else {
		sql = "SELECT ID,boff,blen FROM seq WHERE chrom=? LIMIT 1;";

		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
			sqlite3_bind_text(stmt, 1, name, -1, NULL);
			ret = sqlite3_step(stmt);
		);
		
		if (ret == SQLITE_ROW){
			PYFASTX_SQLITE_CALL(
				chrom = sqlite3_column_int(stmt, 0);
				offset = sqlite3_column_int64(stmt, 1);
				bytes = sqlite3_column_int(stmt, 2);
				sqlite3_finalize(stmt);
			);
		} else {
			PyErr_Format(PyExc_NameError, "Sequence %s does not exists", name);
			PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
			return NULL;
		}

		if (strlen(name) >= self->index->cache_name.m) {
			self->index->cache_name.m = strlen(name) + 1;
			self->index->cache_name.s = (char *)realloc(self->index->cache_name.s, self->index->cache_name.m);
		}

		self->index->cache_full = 1;
		self->index->cache_chrom = chrom;
		strcpy(self->index->cache_name.s, name);
		
		pyfastx_index_fill_cache(self->index, offset, bytes);

		seq = self->index->cache_seq.s;
	}

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
		Py_ssize_t i;
		Py_ssize_t j = 0;
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
	return pyfastx_fasta_keys_create(self->index->index_db, self->seq_counts);
}

PyObject *pyfastx_fasta_subscript(pyfastx_Fasta *self, PyObject *item){
	self->index->iterating = 0;

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
		return pyfastx_index_get_seq_by_name(self->index, item);

	} else {
		PyErr_SetString(PyExc_KeyError, "the key must be index number or sequence name");
		return NULL;
	}
}

Py_ssize_t pyfastx_fasta_length(pyfastx_Fasta *self){
	return self->seq_counts;
}

int pyfastx_fasta_contains(pyfastx_Fasta *self, PyObject *key){
	int ret;
	char *name;
	sqlite3_stmt *stmt;

	if (!PyUnicode_CheckExact(key)) {
		return 0;
	}
	
	name = (char *)PyUnicode_AsUTF8(key);

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, "SELECT 1 FROM seq WHERE chrom=? LIMIT 1;", -1, &stmt, NULL);
		sqlite3_bind_text(stmt, 1, name, -1, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);

	return ret == SQLITE_ROW ? 1 : 0;
}

PyObject *pyfastx_fasta_count(pyfastx_Fasta *self, PyObject *args){
	int ret;

	const char *sql;

	Py_ssize_t l;
	Py_ssize_t c;

	sqlite3_stmt *stmt;

	if (!PyArg_ParseTuple(args, "n", &l)) {
		return NULL;
	}

	sql = "SELECT COUNT(*) FROM seq WHERE slen>=?;";

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		sqlite3_bind_int64(stmt, 1, l);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(c = sqlite3_column_int64(stmt, 0));
	} else {
		c = 0;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return Py_BuildValue("n", c);
}

PyObject *pyfastx_fasta_nl(pyfastx_Fasta *self, PyObject *args){
	int ret;
	int p = 50;
	double half_size;

	const char *sql;

	sqlite3_stmt *stmt;

	Py_ssize_t temp_size = 0;
	Py_ssize_t i = 0;
	Py_ssize_t j = 0;

	if (!PyArg_ParseTuple(args, "|i", &p)) {
		return NULL;
	}

	if (p < 0 || p > 100){
		PyErr_SetString(PyExc_ValueError, "the value must between 0 and 100");
		return NULL;
	}

	if (p == 50) {
		sql = "SELECT n50,l50 FROM stat LIMIT 1";
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
			ret = sqlite3_step(stmt);
		);

		if (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(
				j = sqlite3_column_int64(stmt, 0);
				i = sqlite3_column_int64(stmt, 1);
			);
		}

		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		stmt = NULL;
	}

	if (!j) {
		half_size = p/100.0 * self->seq_length;

		sql = "SELECT slen FROM seq ORDER BY slen DESC";
		PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL));
		for(;;){
			PYFASTX_SQLITE_CALL(ret=sqlite3_step(stmt));
			if (ret == SQLITE_ROW) {
				PYFASTX_SQLITE_CALL(j = sqlite3_column_int64(stmt, 0));
				i++;
				temp_size += j;
				if (temp_size >= half_size)
					break;
			} else {
				break;
			}
		}

		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
	}

	if (j) {
		sql = "UPDATE stat SET n50=?, l50=?";
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
			sqlite3_bind_int64(stmt, 1, j);
			sqlite3_bind_int64(stmt, 2, i);
			sqlite3_step(stmt);
			sqlite3_finalize(stmt);
		);

		return Py_BuildValue("nn", j, i);
	} else {
		PyErr_SetString(PyExc_RuntimeError, "can not calculate N50 and L50");
		return NULL;
	}
}

PyObject *pyfastx_fasta_longest(pyfastx_Fasta *self, void* closure){
	int ret;
	const char *sql;
	sqlite3_stmt *stmt;
	Py_ssize_t chrom;

	sql = "SELECT ID,MAX(slen) FROM seq LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			chrom = sqlite3_column_int64(stmt, 0);
			sqlite3_finalize(stmt);
		);

		return pyfastx_index_get_seq_by_id(self->index, chrom);
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
	PyErr_SetString(PyExc_RuntimeError, "could not found longest sequence");
	return NULL;
}

PyObject *pyfastx_fasta_shortest(pyfastx_Fasta *self, void* closure){
	int ret;
	const char *sql;
	sqlite3_stmt *stmt;
	Py_ssize_t chrom;

	sql = "SELECT ID,MIN(slen) FROM seq LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			chrom = sqlite3_column_int64(stmt, 0);
			sqlite3_finalize(stmt);
		);

		return pyfastx_index_get_seq_by_id(self->index, chrom);
	}
	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
	PyErr_SetString(PyExc_RuntimeError, "not found shortest sequence");
	return NULL;
}

PyObject *pyfastx_fasta_mean(pyfastx_Fasta *self, void* closure){
	int ret;
	double len;
	const char *sql;
	sqlite3_stmt *stmt;

	sql = "SELECT avglen FROM stat LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(len = sqlite3_column_double(stmt, 0));
	} else {
		len = 0;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
	stmt = NULL;

	if (!len) {
		sql = "SELECT AVG(slen) FROM seq LIMIT 1";
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
			ret = sqlite3_step(stmt);
		);

		if (ret == SQLITE_ROW) {
			PYFASTX_SQLITE_CALL(len = sqlite3_column_double(stmt, 0));
		}

		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
	}

	if (len) {
		sql = "UPDATE stat SET avglen=?";
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
			sqlite3_bind_double(stmt, 1, len);
			sqlite3_step(stmt);
			sqlite3_finalize(stmt);
		);

		return Py_BuildValue("d", len);
	} else {
		PyErr_SetString(PyExc_RuntimeError, "could not calculate average length");
		return NULL;
	}
}

PyObject *pyfastx_fasta_median(pyfastx_Fasta *self, void* closure){
	int ret;
	double m;
	const char *sql;
	sqlite3_stmt *stmt;

	sql = "SELECT medlen FROM stat LIMIT 1";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(m = sqlite3_column_double(stmt, 0));
	} else {
		m = 0;
	}

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
	stmt = NULL;

	if (!m) {
		if (self->seq_counts % 2 == 0) {
			sql = "SELECT AVG(slen) FROM (SELECT slen FROM seq ORDER BY slen LIMIT ?,2) LIMIT 1";
		} else {
			sql = "SELECT slen FROM seq ORDER BY slen LIMIT ?,1";
		}

		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
			sqlite3_bind_int64(stmt, 1, (self->seq_counts - 1)/2);
			ret = sqlite3_step(stmt);
		);

		if (ret == SQLITE_ROW){
			PYFASTX_SQLITE_CALL(m = sqlite3_column_double(stmt, 0));
		}
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		stmt = NULL;
	}

	if (m) {
		sql = "UPDATE stat SET medlen=?";
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
			sqlite3_bind_double(stmt, 1, m);
			sqlite3_step(stmt);
			sqlite3_finalize(stmt);
		);

		return Py_BuildValue("d", m);
	} else {
		PyErr_SetString(PyExc_RuntimeError, "could not calculate median length");
		return NULL;
	}
}

PyObject *pyfastx_fasta_format(pyfastx_Fasta *self, void* closure) {
	if (self->index->gzip_format) {
		Py_RETURN_TRUE;
	}

	Py_RETURN_FALSE;
}

void pyfastx_fasta_calc_composition(pyfastx_Fasta *self) {
	int c;
	int j;
	int ret;
	const char *sql;
	sqlite3_stmt *stmt;

	//reading file for kseq
	kstream_t* ks;
	
	//read for line
	kstring_t line = {0, 0, 0};

	//ascii char statistics
	Py_ssize_t seq_comp[128] = {0};
	Py_ssize_t fa_comp[128] = {0};

	Py_ssize_t i;
	Py_ssize_t seqid = 0;

	sql = "SELECT * FROM comp LIMIT 1";

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);
	
	if (ret == SQLITE_ROW)
		return;

	stmt = NULL;

	sql = "PRAGMA synchronous=OFF;BEGIN TRANSACTION;";
	PYFASTX_SQLITE_CALL(sqlite3_exec(self->index->index_db, sql, NULL, NULL, NULL));

	sql = "INSERT INTO comp VALUES (?,?,?,?);";
	PYFASTX_SQLITE_CALL(sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL));
	
	gzrewind(self->index->gzfd);
	ks = ks_init(self->index->gzfd);
	
	Py_BEGIN_ALLOW_THREADS
	
	while (ks_getuntil(ks, '\n', &line, 0) >= 0) {
		if (line.s[0] == 62) {
			if (seqid > 0) {
				for (j = 0; j < 128; ++j) {
					if (seq_comp[j] > 0) {
						sqlite3_bind_null(stmt, 1);
						sqlite3_bind_int64(stmt, 2, seqid);
						sqlite3_bind_int(stmt, 3, j);
						sqlite3_bind_int64(stmt, 4, seq_comp[j]);
						sqlite3_step(stmt);
						sqlite3_reset(stmt);
						fa_comp[j] += seq_comp[j];
					}
				}
			}

			memset(seq_comp, 0, sizeof(seq_comp));
			seqid++;
			continue;
		}

		for (i = 0; i < line.l; ++i) {
			c = line.s[i];

			++seq_comp[c];
		}
	}

	//write the last sequence
	for (j = 0; j < 128; ++j) {
		if (seq_comp[j] > 0) {
			sqlite3_bind_null(stmt, 1);
			sqlite3_bind_int64(stmt, 2, seqid);
			sqlite3_bind_int(stmt, 3, j);
			sqlite3_bind_int64(stmt, 4, seq_comp[j]);
			sqlite3_step(stmt);
			sqlite3_reset(stmt);
			fa_comp[j] += seq_comp[j];
		}
	}

	//write total composition to db
	for (j = 0; j < 128; ++j) {
		sqlite3_bind_null(stmt, 1);
		sqlite3_bind_int64(stmt, 2, 0);
		sqlite3_bind_int(stmt, 3, j);
		sqlite3_bind_int64(stmt, 4, fa_comp[j]);
		sqlite3_step(stmt);
		sqlite3_reset(stmt);
	}

	sqlite3_finalize(stmt);
	sqlite3_exec(self->index->index_db, "COMMIT;", NULL, NULL, NULL);

	Py_END_ALLOW_THREADS

	ks_destroy(ks);
	free(line.s);
}

PyObject *pyfastx_fasta_gc_content(pyfastx_Fasta *self, void* closure) {
	int ret;
	const char *sql;
	sqlite3_stmt *stmt;

	int l;
	Py_ssize_t n;
	Py_ssize_t a = 0;
	Py_ssize_t c = 0;
	Py_ssize_t g = 0;
	Py_ssize_t t = 0;

	pyfastx_fasta_calc_composition(self);
	sql = "SELECT * FROM comp WHERE seqid=0";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

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

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	if (a + c + g + t > 0) {
		return Py_BuildValue("f", (float)(g+c)/(a+c+g+t)*100);
	} else {
		PyErr_SetString(PyExc_RuntimeError, "could not calculate gc content");
		return NULL;
	}
}

PyObject *pyfastx_fasta_gc_skew(pyfastx_Fasta *self, void* closure) {
	int ret;
	const char *sql;
	sqlite3_stmt *stmt;

	int l;
	Py_ssize_t n;
	Py_ssize_t c = 0;
	Py_ssize_t g = 0;

	pyfastx_fasta_calc_composition(self);
	sql = "SELECT * FROM comp WHERE seqid=0";
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

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

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	if (c + g > 0) {
		return Py_BuildValue("f", (float)(g-c)/(g+c));
	} else {
		PyErr_SetString(PyExc_RuntimeError, "could not calculate gc skew");
		return NULL;
	}
}

PyObject *pyfastx_fasta_composition(pyfastx_Fasta *self, void* closure) {
	int l;
	int ret;
	const char *sql;

	sqlite3_stmt *stmt;
	
	Py_ssize_t n;

	PyObject *d;
	PyObject *b;
	PyObject *c;

	pyfastx_fasta_calc_composition(self);

	//the last row store the sum of the each base
	sql = "SELECT * FROM comp WHERE seqid=0";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	d = PyDict_New();

	while (ret == SQLITE_ROW) {
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

	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	return d;
}

//support for guess sequence type according to IUPAC codes
//https://www.bioinformatics.org/sms/iupac.html
PyObject *pyfastx_fasta_guess_type(pyfastx_Fasta *self, void* closure) {
	int l;
	int ret;
	int i;

	char *alphabets;
	char *retval;
	const char *sql;

	sqlite3_stmt *stmt;

	Py_ssize_t n;

	pyfastx_fasta_calc_composition(self);

	sql = "SELECT * FROM comp WHERE seqid=0";
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(self->index->index_db, sql, -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	i = 0;
	alphabets = (char *)malloc(128);

	while (ret == SQLITE_ROW) {
		PYFASTX_SQLITE_CALL(
			l = sqlite3_column_int(stmt, 2);
			n = sqlite3_column_int64(stmt, 3);
			ret = sqlite3_step(stmt);
		);

		if (l > 32 && l < 127 && n > 0) {
			alphabets[i++] = l;
		}
	}

	alphabets[i] = '\0';
	PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));

	if (is_subset("ACGTNacgtn", alphabets) || is_subset("abcdghkmnrstvwyABCDGHKMNRSTVWY*-", alphabets)) {
		retval = "DNA";
	} else if (is_subset("ACGUNacgun", alphabets) || is_subset("abcdghkmnrsuvwyABCDGHKMNRSUVWY*-", alphabets)) {
		retval = "RNA";
	} else if (is_subset("acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY*-", alphabets)) {
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
	{"file_name", T_OBJECT, offsetof(pyfastx_Fasta, file_obj), READONLY},
	{"size", T_PYSSIZET, offsetof(pyfastx_Fasta, seq_length), READONLY},
	{NULL}
};

static PyMethodDef pyfastx_fasta_methods[] = {
	{"build_index", (PyCFunction)pyfastx_fasta_build_index, METH_NOARGS, NULL},
	{"fetch", (PyCFunction)pyfastx_fasta_fetch, METH_VARARGS|METH_KEYWORDS, NULL},
	{"flank", (PyCFunction)pyfastx_fasta_flank, METH_VARARGS|METH_KEYWORDS, NULL},
	{"count", (PyCFunction)pyfastx_fasta_count, METH_VARARGS, NULL},
	{"keys", (PyCFunction)pyfastx_fasta_keys, METH_NOARGS, NULL},
	{"nl", (PyCFunction)pyfastx_fasta_nl, METH_VARARGS, NULL},
	{NULL, NULL, 0, NULL}
};

//as a list
static PySequenceMethods seq_methods = {
	.sq_contains = (objobjproc)pyfastx_fasta_contains,
};

static PyMappingMethods pyfastx_fasta_as_mapping = {
	.mp_length = (lenfunc)pyfastx_fasta_length,
	.mp_subscript = (binaryfunc)pyfastx_fasta_subscript,
};

PyTypeObject pyfastx_FastaType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "Fasta",
	.tp_basicsize = sizeof(pyfastx_Fasta),
	.tp_dealloc = (destructor)pyfastx_fasta_dealloc,
	.tp_repr = (reprfunc)pyfastx_fasta_repr,
	.tp_as_sequence = &seq_methods,
	.tp_as_mapping = &pyfastx_fasta_as_mapping,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_iter = (getiterfunc)pyfastx_fasta_iter,
	.tp_iternext = (iternextfunc)pyfastx_fasta_next,
	.tp_methods = pyfastx_fasta_methods,
	.tp_members = pyfastx_fasta_members,
	.tp_getset = pyfastx_fasta_getsets,
	.tp_new = pyfastx_fasta_new,
};
