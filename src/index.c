#include "index.h"
#include "util.h"
#include "kseq.h"
#include "sequence.h"

/*
create a index
@param file_path, fasta path and name
@param uppercase, uppercase sequence
@param uppercase
*/
pyfastx_Index* pyfastx_init_index(char* file_name, uint16_t uppercase, PyObject* key_func){
	pyfastx_Index* index;

	index = (pyfastx_Index *)malloc(sizeof(pyfastx_Index));
	index->uppercase = uppercase;

	//key function
	index->key_func = key_func;

	//check input file is gzip or not
	index->gzip_format = is_gzip_format(file_name);

	//initial kseqs
	index->gzfd = gzopen(file_name, "rb");
	index->kseqs = kseq_init(index->gzfd);

	//create index file
	index->index_file = (char *)malloc(strlen(file_name) + 5);
	strcpy(index->index_file, file_name);
	strcat(index->index_file, ".fxi");

	index->fd = fopen(file_name, "rb");

	index->index_db = NULL;

	if(index->gzip_format){
		index->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
		//initial zran index
		zran_init(index->gzip_index, index->fd, 4194304, 32768, 1048576, ZRAN_AUTO_BUILD);
	}

	//cache name
	index->cache_chrom = 0;

	//cache start and end position
	index->cache_start = 0;
	index->cache_end = 0;

	//cache sequence
	index->cache_seq = NULL;

	return index;
}

void pyfastx_rewind_index(pyfastx_Index *index){
	kseq_rewind(index->kseqs);
	gzrewind(index->gzfd);
}

PyObject* pyfastx_get_next_seq(pyfastx_Index *index){
	if(kseq_read(index->kseqs) >= 0){
		if(index->uppercase){
			upper_string(index->kseqs->seq.s);
		}
		return Py_BuildValue("(ss)", index->kseqs->name.s, index->kseqs->seq.s);
	}
	return NULL;
}

void pyfastx_create_index(pyfastx_Index *self){
	// seqlite3 return value
	//int ret;
	
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	
	// 1: normal fasta sequence with the same length in line
	// 0: not normal fasta sequence with different length in line
	uint16_t seq_normal = 1;

	//1: \n, 2: \r\n
	uint16_t line_end = 1;
	
	//length of previous line
	uint32_t line_len = 0;

	//length of current line
	uint32_t temp_len = 0;

	//number of lines that line_len not equal to temp_len
	uint32_t bad_line = 0;

	//current read position
	uint64_t position = 0;
	
	// start position
	uint64_t start = 0;

	//sequence length
	uint32_t seq_len = 0;

	//number of bases
	uint32_t g_count = 0;
	uint32_t c_count = 0;
	uint32_t a_count = 0;
	uint32_t t_count = 0;
	uint32_t n_count = 0;

	//current read base char
	uint32_t i;

	//real line len
	uint32_t real_len;

	//reading file for kseq
	kstream_t* ks;

	//read for line
	kstring_t line = {0, 0, 0};

	//description line
	char* description = NULL;

	//chromosome name
	char* chrom = NULL;

	if (sqlite3_open(self->index_file, &self->index_db) != SQLITE_OK) {
		PyErr_SetString(PyExc_ConnectionError, sqlite3_errmsg(self->index_db));
		return;
	}

	//create index database
	const char *create_sql = " \
		CREATE TABLE seq ( \
			ID INTEGER PRIMARY KEY, --seq identifier\n \
			chrom TEXT, --seq name\n \
			boff INTEGER, --seq offset start\n \
			blen INTEGER, --seq byte length\n \
			slen INTEGER, --seq length\n \
			llen INTEGER, --line lenght\n \
			elen INTEGER, --end length\n \
			norm INTEGER, --line with the same length or not\n \
			a INTEGER, --A base counts\n \
			c INTEGER, --C base counts\n \
			g INTEGER, --G base counts\n \
			t INTEGER, --T base counts\n \
			n INTEGER, --unknown base counts\n \
			descr TEXT --sequence description\n \
		); \
		CREATE TABLE gzindex ( \
			ID INTEGER PRIMARY KEY, \
			content BLOB \
		);";

	if(sqlite3_exec(self->index_db, create_sql, NULL, NULL, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return;
	}

	if(sqlite3_exec(self->index_db, "PRAGMA synchronous=OFF;BEGIN;", NULL, NULL, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return;
	}

	const char *insert_sql = "INSERT INTO seq VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	sqlite3_prepare_v2(self->index_db, insert_sql, -1, &stmt, NULL);

	Py_BEGIN_ALLOW_THREADS

	ks = ks_init(self->gzfd);

	while (ks_getuntil(ks, '\n', &line, 0) >= 0) {
		position += line.l;
		position++;

		if (line.s[0] == '>') {
			if (start > 0) {
				//end of sequence and check whether normal fasta
				seq_normal = (bad_line > 1) ? 0 : 1;

				sqlite3_bind_null(stmt, 1);
				sqlite3_bind_text(stmt, 2, chrom, -1, NULL);
				sqlite3_bind_int64(stmt, 3, start);
				sqlite3_bind_int(stmt, 4, position-start-line.l-1);
				sqlite3_bind_int(stmt, 5, seq_len);
				sqlite3_bind_int(stmt, 6, line_len);
				sqlite3_bind_int(stmt, 7, line_end);
				sqlite3_bind_int(stmt, 8, seq_normal);
				sqlite3_bind_int(stmt, 9, a_count);
				sqlite3_bind_int(stmt, 10, c_count);
				sqlite3_bind_int(stmt, 11, g_count);
				sqlite3_bind_int(stmt, 12, t_count);
				sqlite3_bind_int(stmt, 13, n_count);
				sqlite3_bind_text(stmt, 14, description, -1, NULL);
				sqlite3_step(stmt);
				sqlite3_reset(stmt);
			}

			//reset
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
			bad_line = 0;
			seq_normal = 1;

			//get line end length \r\n or \n
			if (line.s[line.l-1] == '\r') {
				line_end = 2;
			}

			description = (char *)malloc(line.l);
			memcpy(description, line.s+1, line.l-line_end);
			description[line.l-line_end] = '\0';
			

			if (self->key_func == Py_None) {
				chrom = (char *)malloc(line.l);
				strcpy(chrom, description);
				strtok(chrom, " ");
			} else {
				PyGILState_STATE state = PyGILState_Ensure();
				PyObject *result = PyObject_CallFunction(self->key_func, "s", description);
				PyGILState_Release(state);
				chrom = PyUnicode_AsUTF8(result);
			}

			continue;
		}

		temp_len = line.l + 1;

		if (line_len > 0 && line_len != temp_len) {
			bad_line++;
		}

		//record first line length
		if (line_len == 0) {
			line_len = temp_len;
		}

		//calculate atgc counts
		real_len = line.l - line_end + 1;

		//calculate seq len
		seq_len += real_len;

		for (i = 0; i < real_len; i++) {
			switch (line.s[i]) {
				case 65: case 97: ++a_count; break;
				case 67: case 99: ++c_count; break;
				case 71: case 103: ++g_count; break;
				case 84: case 116: ++t_count; break;
				default: ++n_count; break;
			}
		}
	}

	//end of sequence and check whether normal fasta
	seq_normal = (bad_line > 1) ? 0 : 1;

	sqlite3_bind_null(stmt, 1);
	sqlite3_bind_text(stmt, 2, chrom, -1, NULL);
	sqlite3_bind_int64(stmt, 3, start);
	sqlite3_bind_int(stmt, 4, position-start);
	sqlite3_bind_int(stmt, 5, seq_len);
	sqlite3_bind_int(stmt, 6, line_len);
	sqlite3_bind_int(stmt, 7, line_end);
	sqlite3_bind_int(stmt, 8, seq_normal);
	sqlite3_bind_int(stmt, 9, a_count);
	sqlite3_bind_int(stmt, 10, c_count);
	sqlite3_bind_int(stmt, 11, g_count);
	sqlite3_bind_int(stmt, 12, t_count);
	sqlite3_bind_int(stmt, 13, n_count);
	sqlite3_bind_text(stmt, 14, description, -1, NULL);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);

	sqlite3_exec(self->index_db, "CREATE INDEX chromidx ON seq (chrom);", NULL, NULL, NULL);
	sqlite3_exec(self->index_db, "COMMIT;", NULL, NULL, NULL);

	ks_destroy(ks);
	free(line.s);

	//create gzip random access index
	if(self->gzip_format){
		pyfastx_build_gzip_index(self->gzip_index, self->index_db, self->index_file);
	}

	Py_END_ALLOW_THREADS
}

//load index from index file
void pyfastx_load_index(pyfastx_Index *self){
	if(sqlite3_open(self->index_file, &self->index_db) != SQLITE_OK){
		PyErr_SetString(PyExc_ConnectionError, sqlite3_errmsg(self->index_db));
		return;
	}

	if(self->gzip_format){
		pyfastx_load_gzip_index(self->gzip_index, self->index_db, self->index_file);
	}
}

void pyfastx_build_index(pyfastx_Index *self){
	if(file_exists(self->index_file)) {
		pyfastx_load_index(self);
	} else {
		pyfastx_create_index(self);
	}
}

void pyfastx_index_free(pyfastx_Index *self){
	if (self->gzip_format) {
		zran_free(self->gzip_index);
	}

	if (self->index_db != NULL) {
		sqlite3_close(self->index_db);
	}

	if (self->cache_seq != NULL) {
		free(self->cache_seq);
	}

	kseq_destroy(self->kseqs);
	gzclose(self->gzfd);
	fclose(self->fd);
}

PyObject *pyfastx_index_make_seq(pyfastx_Index *self, sqlite3_stmt *stmt){
	int32_t a, c, g, t, n;
	char* name;

	pyfastx_Sequence *seq = PyObject_New(pyfastx_Sequence, &pyfastx_SequenceType);
	if(!seq){
		return NULL;
	}

	seq->id = sqlite3_column_int(stmt, 0);
	name = (char *)sqlite3_column_text(stmt, 1);
	seq->name = (char *)malloc(strlen(name) + 1);
	strcpy(seq->name, name);
	seq->offset = (int64_t)sqlite3_column_int64(stmt, 2);
	seq->byte_len = sqlite3_column_int(stmt, 3);
	seq->seq_len = sqlite3_column_int(stmt, 4);
	seq->parent_len = seq->seq_len;
	seq->line_len = sqlite3_column_int(stmt, 5);
	seq->end_len = sqlite3_column_int(stmt, 6);
	seq->normal = sqlite3_column_int(stmt, 7);
	a = sqlite3_column_int(stmt, 8);
	c = sqlite3_column_int(stmt, 9);
	g = sqlite3_column_int(stmt, 10);
	t = sqlite3_column_int(stmt, 11);
	n = sqlite3_column_int(stmt, 12);

	sqlite3_finalize(stmt);

	//position
	seq->start = 1;
	seq->end = seq->seq_len;

	//composition
	seq->composition = Py_BuildValue("{s:I,s:I,s:I,s:I,s:I}", "A", a, "C", c, "G", g, "T", t, "N", n);

	//calc GC content
	seq->gc_content = (float)(g+c)/(a+c+g+t)*100;

	//calc GC skew
	seq->gc_skew = (float)(g-c)/(g+c);

	//index
	seq->index = self;

	Py_INCREF(seq);
	return (PyObject *)seq;
}

PyObject *pyfastx_index_get_seq_by_name(pyfastx_Index *self, char *name){
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	
	//select sql statement, chrom indicates seq name or chromomsome
	const char* sql = "SELECT * FROM seq WHERE chrom=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, name, -1, NULL);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		PyErr_SetString(PyExc_KeyError, name);
		return NULL;
	}

	return pyfastx_index_make_seq(self, stmt);
}


PyObject *pyfastx_index_get_seq_by_id(pyfastx_Index *self, uint32_t chrom){
	sqlite3_stmt *stmt;

	const char* sql = "SELECT * FROM seq WHERE ID=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, chrom);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		PyErr_SetString(PyExc_IndexError, "Index Error");
		return NULL;
	}

	return pyfastx_index_make_seq(self, stmt);
}

char *pyfastx_index_get_full_seq(pyfastx_Index *self, uint32_t chrom){
	sqlite3_stmt *stmt;
	uint32_t seq_len;
	int64_t offset;
	uint32_t bytes;
	char *buff;
	
	//select sql statement, chrom indicates seq or chromomsome id
	const char* sql = "SELECT boff,blen,slen FROM seq WHERE ID=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, chrom);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		PyErr_SetString(PyExc_KeyError, "Can not found sequence");
		return NULL;
	}

	offset = (int64_t)sqlite3_column_int64(stmt, 0);
	bytes = sqlite3_column_int(stmt, 1);
	seq_len = sqlite3_column_int(stmt, 2);

	if((chrom == self->cache_chrom) && (1==self->cache_start) && (seq_len==self->cache_end)){
		return self->cache_seq;
	}

	Py_BEGIN_ALLOW_THREADS
	
	if(self->gzip_format){
		buff = (char *)malloc(bytes + 1);
		zran_seek(self->gzip_index, offset, SEEK_SET, NULL);
		zran_read(self->gzip_index, buff, bytes);
		buff[bytes] = '\0';
		remove_space(buff);
	} else {
		gzseek(self->gzfd, offset, SEEK_SET);
		kstream_t *ks;
		kstring_t seq = {0, 0, 0};
		seq.m = 256;
		seq.s = (char*)malloc(seq.m);

		ks = ks_init(self->gzfd);
		
		int c;
		while((c = ks_getc(ks)) >= 0 && c != '>'){
			if(c == '\n') continue;
			seq.s[seq.l++] = c;
			ks_getuntil2(ks, 2, &seq, 0, 1);
		}

		seq.s[seq.l] = 0;
		buff = seq.s;
	}
	
	if(self->uppercase) {
		upper_string(buff);
	}

	Py_END_ALLOW_THREADS

	self->cache_chrom = chrom;
	self->cache_start = 1;
	self->cache_end = seq_len;
	self->cache_seq = buff;

	return self->cache_seq;
}

/*
@name str, sequence identifier in fasta file
@offset int, subsequence byte start in fasta file
@bytes int, byte length of subsequence with space
@start int, sliced start position in sequence
@end int, sliced end position in sequence
@plen int, length of subsequence parent seq
@normal int, 1 -> normal fasta with the same length, 0 not normal
@return str
*/
char *pyfastx_index_get_sub_seq(pyfastx_Index *self, uint32_t chrom, int64_t offset, int64_t bytes, uint32_t start, uint32_t end, uint32_t plen, uint16_t normal){
	uint32_t seq_len;
	char *buff;
	seq_len = end - start + 1;

	if(!normal || (plen == end && start == 1)) {
		buff = pyfastx_index_get_full_seq(self, chrom);
	}
	
	if((chrom == self->cache_chrom) && (start==self->cache_start) && (end==self->cache_end)){
		return self->cache_seq;
	}

	if((chrom == self->cache_chrom) && (start>=self->cache_start) && (end<=self->cache_end)){
		buff = (char *)malloc(seq_len + 1);
		memcpy(buff, self->cache_seq + (start - self->cache_start), seq_len);
		buff[seq_len] = '\0';
		return buff;
	}
	
	buff = (char *)malloc(bytes + 1);

	Py_BEGIN_ALLOW_THREADS

	if(self->gzip_format){
		zran_seek(self->gzip_index, offset, SEEK_SET, NULL);
		zran_read(self->gzip_index, buff, bytes);
	} else {
		gzseek(self->gzfd, offset, SEEK_SET);
		gzread(self->gzfd, buff, bytes);
	}

	buff[bytes] = '\0';

	remove_space(buff);

	if(self->uppercase){
		upper_string(buff);
	}

	Py_END_ALLOW_THREADS

	self->cache_chrom = chrom;
	self->cache_start = start;
	self->cache_end = end;
	self->cache_seq = buff;


	return buff;
}
