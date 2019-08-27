#include "index.h"
#include "util.h"
#include "sequence.h"

/*
create a index
@param file_path, fasta path and name
@param uppercase, uppercase sequence
@param uppercase
*/
pyfastx_Index* pyfastx_init_index(char* file_name, int uppercase){
	pyfastx_Index* index;

	index = (pyfastx_Index *)malloc(sizeof(pyfastx_Index));
	index->uppercase = uppercase;

	//check input file is gzip or not
	index->gzip_format = is_gzip_format(file_name);

	//initial kseqs
	index->gzfd = gzopen(file_name, "rb");
	index->kseqs = kseq_init(index->gzfd);

	//create index file
	index->index_file = (char *)malloc(strlen(file_name) + 4);
	strcpy(index->index_file, file_name);
	strcat(index->index_file, ".db");

	index->fd = fopen(file_name, "rb");

	index->index_db = NULL;

	if(index->gzip_format){
		index->gzip_index = (zran_index_t *)malloc(sizeof(zran_index_t));
	}

	//cache name
	index->cache_name = NULL;

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

void pyfastx_build_gzip_index(pyfastx_Index *self){
	sqlite3_stmt *stmt;

	rewind(self->fd);
	zran_init(self->gzip_index, self->fd, 4194304, 32768, 1048576, ZRAN_AUTO_BUILD);
	zran_build_index(self->gzip_index, 0, 0);

	//create temp gzip index file
	char *temp_index = (char *)malloc(strlen(self->index_file) + 5);
	strcpy(temp_index, self->index_file);
	strcat(temp_index, ".tmp");
	
	FILE* fh = fopen(temp_index, "wb+");
	FILE* fd = fdopen(fileno(fh), "ab");
	
	zran_export_index(self->gzip_index, fd);
	
	int fsize = ftell(fh);
	rewind(fh);

	char *buff = (char *)malloc(fsize + 1);

	if(fread(buff, fsize, 1, fh) != 1){
		return;
	}
	
	fclose(fh);
	remove(temp_index);

	sqlite3_prepare_v2(self->index_db, "INSERT INTO gzindex VALUES (NULL, ?)", -1, &stmt, NULL);
	sqlite3_bind_blob(stmt, 1, buff, fsize, NULL);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
	free(buff);
}

void pyfastx_load_gzip_index(pyfastx_Index *self){
	sqlite3_stmt *stmt;
	int bytes = 0;
	
	zran_init(self->gzip_index, self->fd, 4194304, 32768, 1048576, ZRAN_AUTO_BUILD);
	
	char *temp_index = (char *)malloc(strlen(self->index_file) + 5);
	strcpy(temp_index, self->index_file);
	strcat(temp_index, ".tmp");
	
	FILE* fh = fopen(temp_index, "wb+");

	sqlite3_prepare_v2(self->index_db, "SELECT content FROM gzindex;", -1, &stmt, NULL);
	if(sqlite3_step(stmt) == SQLITE_ROW){
		bytes = sqlite3_column_bytes(stmt, 0);
	}
	//const char *buff = sqlite3_column_blob(stmt, 0);
	//fwrite(buff, 1, strlen(buff), fh);
	fwrite(sqlite3_column_blob(stmt, 0), bytes, 1, fh);
	//fseek(fh, 0, SEEK_SET);
	rewind(fh);

	FILE* fd = fdopen(fileno(fh), "rb");
	zran_import_index(self->gzip_index, fd);

	fclose(fh);
	remove(temp_index);
}

void pyfastx_create_index(pyfastx_Index *self){
	// seqlite3 return value
	//int ret;
	
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
	kstream_t* ks = self->kseqs->f;

	if(sqlite3_open(self->index_file, &self->index_db) != SQLITE_OK){
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
			norm INTEGER, --line with the same length or not\n \
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

	if(sqlite3_exec(self->index_db, create_sql, NULL, NULL, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return;
	}

	if(sqlite3_exec(self->index_db, "PRAGMA synchronous=OFF;BEGIN;", NULL, NULL, NULL) != SQLITE_OK){
		PyErr_SetString(PyExc_RuntimeError, sqlite3_errmsg(self->index_db));
		return;
	}

	const char *insert_sql = "INSERT INTO seq VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?);";
	sqlite3_prepare_v2(self->index_db, insert_sql, -1, &stmt, NULL);

	
	while((c=ks_getc(ks)) >= 0){
		position++;
	
		switch(c){
			// c is >
			case 62: {
				if(start > 0){

					//end of sequence and check whether normal fasta
					seq_normal = (bad_line > 1) ? 0 : 1;

					sqlite3_bind_null(stmt, 1);
					sqlite3_bind_text(stmt, 2, self->kseqs->name.s, self->kseqs->name.l, NULL);
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
				bad_line = 0;
				seq_normal = 1;

				break;
			}

			// c is \r
			case 13: {
				temp_len++;
				if(line_end != 2){
					line_end = 2;
				}
				break;
			}
			
			// c is \n
			case 10: {
				temp_len++;
				
				if(line_len > 0 && line_len != temp_len){
					bad_line++;
				}
				if(line_len == 0){
					line_len = temp_len;
				}
				temp_len = 0;
				break;
			}

			default: {
				seq_len++;

				//temp line length
				temp_len++;

				//c = toupper(c);
				
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

					case 97:
						a_count++; break;

					case 99:
						c_count++; break;

					case 103:
						g_count++; break;

					case 116:
						t_count++; break;

					default:
						n_count++;
				}
			}
		}
	}

	//end of sequenc and check whether normal fasta
	seq_normal = (bad_line > 1) ? 0 : seq_normal;

	sqlite3_bind_null(stmt, 1);
	sqlite3_bind_text(stmt, 2, self->kseqs->name.s, self->kseqs->name.l, NULL);
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
		pyfastx_build_gzip_index(self);
	}
}

//load index from index file
void pyfastx_load_index(pyfastx_Index *self){
	if(sqlite3_open(self->index_file, &self->index_db) != SQLITE_OK){
		PyErr_SetString(PyExc_ConnectionError, sqlite3_errmsg(self->index_db));
		return;
	}

	if(self->gzip_format){
		pyfastx_load_gzip_index(self);
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
	if(self->gzip_format){
		zran_free(self->gzip_index);
	}
	if(self->index_db != NULL){
		sqlite3_close(self->index_db);
	}
	kseq_destroy(self->kseqs);
	gzclose(self->gzfd);
}

PyObject *pyfastx_index_make_seq(pyfastx_Index *self, sqlite3_stmt *stmt){
	int a, c, g, t, n;
	char* name;

	pyfastx_Sequence *seq = PyObject_New(pyfastx_Sequence, &pyfastx_SequenceType);
	if(!seq){
		return NULL;
	}

	//seq->index_id = sqlite3_column_int(stmt, 0);
	name = (char *) sqlite3_column_text(stmt, 1);
	seq->name = (char *)malloc(strlen(name) + 1);
	strcpy(seq->name, name);
	seq->offset = (int64_t)sqlite3_column_int64(stmt, 2);
	seq->byte_len = sqlite3_column_int(stmt, 3);
	seq->seq_len = sqlite3_column_int(stmt, 4);
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
	seq->composition = Py_BuildValue("{s:i,s:i,s:i,s:i,s:i}", "A", a, "C", c, "G", g, "T", t, "N", n);

	//calc GC content
	seq->gc_content = (float)(g+c)/(a+c+g+t)*100;

	//index
	seq->index = self;

	Py_INCREF(seq);
	return (PyObject *)seq;
}

PyObject *pyfastx_index_get_seq_by_name(pyfastx_Index *self, char *name){
	// sqlite3 prepare object
	sqlite3_stmt *stmt;
	
	//select sql statement, seqid indicates seq name or chromomsome
	const char* sql = "SELECT * FROM seq WHERE seqid=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, name, -1, NULL);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		PyErr_SetString(PyExc_KeyError, name);
		return NULL;
	}

	return pyfastx_index_make_seq(self, stmt);
}


PyObject *pyfastx_index_get_seq_by_id(pyfastx_Index *self, int id){
	sqlite3_stmt *stmt;

	const char* sql = "SELECT * FROM seq WHERE id=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_int(stmt, 1, id);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		PyErr_SetString(PyExc_IndexError, "Index Error");
		return NULL;
	}

	return pyfastx_index_make_seq(self, stmt);
}

char *pyfastx_index_get_full_seq(pyfastx_Index *self, char *name){
	sqlite3_stmt *stmt;
	int seq_len;
	int64_t offset;
	int bytes;
	
	//select sql statement, seqid indicates seq name or chromomsome
	const char* sql = "SELECT * FROM seq WHERE seqid=? LIMIT 1;";
	sqlite3_prepare_v2(self->index_db, sql, -1, &stmt, NULL);
	sqlite3_bind_text(stmt, 1, name, -1, NULL);
	if(sqlite3_step(stmt) != SQLITE_ROW){
		PyErr_SetString(PyExc_KeyError, name);
		return NULL;
	}

	offset = (int64_t)sqlite3_column_int64(stmt, 2);
	bytes = sqlite3_column_int(stmt, 3);
	seq_len = sqlite3_column_int(stmt, 4);

	if((name==self->cache_name) && (1==self->cache_start) && (seq_len==self->cache_end)){
		return self->cache_seq;
	}

	char *buff = (char *)malloc(bytes + 1);
	if(self->gzip_format){
		zran_seek(self->gzip_index, offset, SEEK_SET, NULL);
		zran_read(self->gzip_index, buff, bytes);
	} else {
		fseek(self->fd, offset, SEEK_SET);
		if(fread(buff, bytes, 1, self->fd) != 1){
			return NULL;
		}
	}

	buff[bytes] = '\0';

	remove_space(buff);

	if(self->uppercase){
		upper_string(buff);
	}

	self->cache_name = name;
	self->cache_start = 1;
	self->cache_end = seq_len;
	self->cache_seq = buff;

	return self->cache_seq;
}

/*
@name str, sequence identifier in fasta file
@offset int, sequence byte start in fasta file
@bytes int, byte length of sequence with space
@start int, start position in sequence
@end int, end position in sequence
@normal int, 1 -> normal fasta with the same length, 0 not normal
@return str
*/
char *pyfastx_index_get_sub_seq(pyfastx_Index *self, char *name, int64_t offset, int bytes, int start, int end, int normal){
	int seq_len;
	char *buff;
	seq_len = end - start + 1;

	if(!normal) {
		buff = pyfastx_index_get_full_seq(self, name);
	}

	if(self->cache_name != NULL){
		if((strcmp(name,self->cache_name)==0) && (start==self->cache_start) && (end==self->cache_end)){
			return self->cache_seq;
		}

		if((strcmp(name,self->cache_name)==0) && (start>=self->cache_start) && (end<=self->cache_end)){
			buff = (char *)malloc(seq_len + 1);
			strncpy(buff, self->cache_seq + (start - self->cache_start), seq_len);
			buff[seq_len] = '\0';
			return buff;
		}
	}


	buff = (char *)malloc(bytes + 1);

	if(self->gzip_format){
		zran_seek(self->gzip_index, offset, SEEK_SET, NULL);
		zran_read(self->gzip_index, buff, bytes);
	} else {
		//rewind(self->fd);
		fseek(self->fd, offset, SEEK_SET);
		if(fread(buff, bytes, 1, self->fd) != 1){
			return NULL;
		}
	}

	buff[bytes] = '\0';

	remove_space(buff);

	if(self->uppercase){
		upper_string(buff);
	}

	self->cache_name = name;
	self->cache_start = start;
	self->cache_end = end;
	self->cache_seq = buff;

	return buff;
}
