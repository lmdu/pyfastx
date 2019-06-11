#include "index.h"
#include "util.h"
//#include "structmember.h"
KSEQ_INIT(gzFile, gzread, gzrewind)

static void _pyfastx_build_gzip_index(pyfastx_Fasta *self){
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

	sqlite3_prepare_v2(self->index_db, "INSERT INTO gzindex VALUES (NULL, ?)", -1, &stmt, NULL);
	sqlite3_bind_blob(stmt, 1, buff, strlen(buff), NULL);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
}

static void _pyfastx_load_gzip_index(pyfastx_Fasta *self){
	sqlite3_stmt *stmt;
	zran_init(self->gzip_index, self->fd, 0, 0, 0, ZRAN_AUTO_BUILD);
	char *temp_index = malloc(strlen(self->file_name) + 5);
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
static void _pyfastx_calc_fasta_attrs(pyfastx_Fasta *self){
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
	sqlite3_prepare_v2(self->index_db, "SELECT SUM(slen) FROM seq LIMIT 1;", -1 &stmt, NULL);
	sqlite3_step(stmt);
	self->seq_length = sqlite3_column_int64(stmt, 0);
	sqlite3_reset(stmt);

	//calculate base counts
	sqlite3_prepare_v2(self->index_db, "SELECT SUM(a),SUM(c),SUM(g),SUM(t),SUM(n) FROM seq LIMIT 1;", -1 &stmt, NULL);
	sqlite3_step(stmt);
	a_counts = sqlite3_column_int(stmt, 0);
	c_counts = sqlite3_column_int(stmt, 1);
	g_counts = sqlite3_column_int(stmt, 2);
	t_counts = sqlite3_column_int(stmt, 3);
	n_counts = sqlite3_column_int(stmt, 4);
	self->composition = Py_BuildValue("{s:i,s:i,s:i,s:i,s:i}", "A", a_counts, "C", c_counts, "G", g_counts, "T", t_counts, "N", n_counts);
	sqlite3_finalize(stmt);

	//calc GC content
	self->gc_content = ((g_counts+c_counts)*1.0/(a_counts+c_counts+g_counts+t_counts)*100*100+0.5)/100.0;
}

PyObject *pyfastx_build_index(pyfastx_Fasta *self, PyObject *args, PyObject *kwargs){
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
		}
	}

	ret = sqlite3_open(self->index_file, &self->index_db);
	if(ret != SQLITE_OK){
		return Py_BuildValue("i", 0);
	}

	//create index database
	const char *create_sql = " \
		CREATE TABLE seq ( \
			sid TEXTPRIMARY KEY, --seq id\n \
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

	ret = sqlite3_exec(self->index_db, "PRAGMA synchronous=OFF", NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "pragma");
	}

	ret = sqlite3_exec(self->index_db, "begin", NULL, NULL, NULL);
	if(ret != SQLITE_OK){
		printf("%s\n", "begin error");
	}

	char *insert_sql = "INSERT INTO seq VALUES (?,?,?,?,?,?,?,?,?,?,?,?)";
	
	sqlite3_prepare_v2(self->index_db, insert_sql, -1, &stmt, NULL);

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

	sqlite3_exec(self->index_db, "commit;", NULL, NULL, NULL);

	//create gzip random access index
	if(self->gzip_format){
		_pyfastx_build_gzip_index(self);
	}

	//get attributes
	_pyfastx_calc_fasta_attrs(self);
	
	return Py_BuildValue("i", 1);
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
