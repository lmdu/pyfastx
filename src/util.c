#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "util.h"

#ifdef _WIN32
#include "fcntl.h"
#include "stdio.h"

int mkstemp(char *template) {
	if (_mktemp_s(template, strlen(template) + 1) != 0) {
		return -1;
	}

	return open(template, _O_CREAT | _O_EXCL, _S_IWRITE);
}
#endif

//check file is whether exists in disk
int file_exists(char *file_name){
	FILE *file;

	if((file = fopen(file_name, "r"))){
		fclose(file);
		return 1;
	}

	return 0;
}

//check file is fasta file
int fasta_validator(gzFile fd) {
	int c;

	while ((c=gzgetc(fd)) != -1) {
		if (isspace(c)) {
			continue;
		}

		if (c == '>') {
			return 1;
		} else {
			return 0;
		}
	}

	return 0;
}

int fastq_validator(gzFile fd) {
	int c;

	while ((c=gzgetc(fd)) != -1) {
		if (isspace(c)) {
			continue;
		}

		if (c == '@') {
			return 1;
		} else {
			return 0;
		}
	}

	return 0;
}

//check file is fasta or fastq file
int fasta_or_fastq(gzFile fd) {
	int c;

	while ((c=gzgetc(fd)) != -1) {
		if (isspace(c)) {
			continue;
		}

		if (c == '>') {
			return 1;
		} else if ( c == '@') {
			return 2;
		} else {
			return 0;
		}
	}

	return 0;
}

/*
remove space using jump table
Reference:
https://github.com/lemire/despacer/blob/master/include/despacer.h
*/
int jump_table[128] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
};

Py_ssize_t remove_space(char *str, Py_ssize_t len) {
	unsigned char c;
	Py_ssize_t i = 0, j = 0;

	while (i < len) {
		c = str[i++];
		str[j] = c;
		j += jump_table[c];
	}

	str[j] = '\0';

	return j;
}

Py_ssize_t remove_space_uppercase(char *str, Py_ssize_t len) {
	unsigned char c;
	Py_ssize_t i = 0, j = 0;

	while (i < len) {
		c = str[i++];
		str[j] = Py_TOUPPER(Py_CHARMASK(c));
		j += jump_table[c];
	}

	str[j] = '\0';

	return j;
}

void upper_string(char *str, Py_ssize_t len) {
	Py_ssize_t i;

	for (i = 0; i < len; ++i) {
		str[i] = Py_TOUPPER(Py_CHARMASK(str[i]));
	}
}

/* 
DNA IUPAC Ambiguity Codes
IUPAC Codes	Meaning		Complement
A			A			T
C			C			G
G			G			C
T/U			T			A
M			A,C			K
R			A,G			Y
W			A,T			W
S			C,G			S
Y			C,T			R
K			G,T			M
V			A,C,G		B
H			A,C,T		D
D			A,G,T		H
B			C,G,T		V
N			G,A,T,C		N

References:
https://droog.gs.washington.edu/parc/images/iupac.html
http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
https://github.com/samtools/samtools/blob/f6bd3b22ea4c9e4897cda455e786180fe650e494/faidx.c
*/
int comp_map[128] = {
	  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
	 16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
	 32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
	 48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
	 64,  84,  86,  71,  72,  69,  70,  67,  68,  73,  74,  77,  76,  75,  78,  79,
	 80,  81,  89,  83,  65,  65,  66,  87,  88,  82,  90,  91,  92,  93,  94,  95,
	 96, 116, 118, 103, 104, 101, 102,  99, 100, 105, 106, 109, 108, 107, 110, 111,
	112, 113, 121, 115,  97,  97,  98, 119, 120, 114, 122, 123, 124, 125, 126, 127,
};

void reverse_complement_seq(char *seq) {
	char c;
	char *p1 = seq;
	char *p2 = seq + strlen(seq) - 1;

	while (p1 <= p2) {
		c = comp_map[Py_CHARMASK(*p1)];
		*p1++ = comp_map[Py_CHARMASK(*p2)];
		*p2-- = c;
	}
}

void reverse_seq(char *seq) {
	int c;
	char *p1 = seq;
	char *p2 = seq + strlen(seq) - 1;

	while (p1 < p2) {
		c = *p1;
		*p1++ = *p2;
		*p2-- = c;
	}
}

void complement_seq(char *seq) {
	int c;

	while ((c = *seq)) {
		*seq++ = comp_map[Py_CHARMASK(c)];
	}
}

Py_ssize_t sum_array(Py_ssize_t arr[], int num) {
	int i;
	Py_ssize_t sum=0;

	for (i = 0; i < num; ++i) {
		sum += arr[i];
	}

	return sum;
}

int is_subset(char *seq1, char *seq2) {
	size_t i, j, m, n;

	m = strlen(seq1);
	n = strlen(seq2);

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (seq2[i] == seq1[j]) {
				break;
			}
		}

		if (j == m) {
			return 0;
		}
	}

	return 1;
}

/*check input file is whether gzip file
@para file_name str, input file path string
@return bool, 1 is gzip formmat file, 0 is not gzip
*/
int is_gzip_format(char* file_name){
	int ret;
	FILE* fd;
	unsigned char magic[4] = {0};

	fd = fopen(file_name, "rb");
	ret = fread(magic, sizeof(magic), 1, fd);
	fclose(fd);

	if (ret != 1){
		return 0;
	}
	
	if (magic[0] != 0x1f || magic[1] != 0x8b || magic[2] != 0x08){
		return 0;
	}
	
	return 1;
}

//generate random file name according to given length
char *generate_random_name(char* index_file) {
	int file_len;
	char *result;

	file_len = strlen(index_file);
	result = (char *)malloc(file_len + 8);
	sprintf(result, "%s.XXXXXX", index_file);

	return result;
}

void pyfastx_build_gzip_index(char* index_file, zran_index_t* gzip_index, sqlite3* index_db) {
	int fd;
	int ret;
	int rowid;

	char *temp_index;
	void *buff;

	FILE* fh;

	sqlite3_stmt *stmt;
	sqlite3_blob *blob;
	
	Py_ssize_t remain;
	Py_ssize_t offset;
	Py_ssize_t block;
	Py_ssize_t len;
	Py_ssize_t chunk;

	ret = zran_build_index(gzip_index, 0, 0);

	if (ret != 0) {
		PyErr_SetString(PyExc_RuntimeError, "failed to build gzip index");
		return;
	}

	temp_index = generate_random_name(index_file);
	if ((fd = mkstemp(temp_index)) < 0) {
		PyErr_SetString(PyExc_RuntimeError, "failed to create temp file");
		return;
	}
	close(fd);

	fh = fopen(temp_index, "wb+");
	if (zran_export_index(gzip_index, fh, NULL) != ZRAN_EXPORT_OK){
		fclose(fh);
		free(temp_index);
		PyErr_SetString(PyExc_RuntimeError, "failed to export gzip index");
		return;
	}

	remain = FTELL(fh);
	rewind(fh);

	buff = malloc(1048576);

	while (remain > 0) {
		if (remain > 524288000) {
			block = 524288000;
		} else {
			block = remain;
		}
		offset = 0;
		PYFASTX_SQLITE_CALL(
			sqlite3_prepare_v2(index_db, "INSERT INTO gzindex VALUES (?,?)", -1, &stmt, NULL);
			sqlite3_bind_null(stmt, 1);
			sqlite3_bind_zeroblob(stmt, 2, block);
			sqlite3_step(stmt);
			rowid = sqlite3_last_insert_rowid(index_db);
			sqlite3_blob_open(index_db, "main", "gzindex", "content", rowid, 1, &blob);

			while (offset < block) {
				chunk = block - offset;
				if (chunk > 1048576) chunk = 1048576;

				if ((len=fread(buff, 1, chunk, fh)) > 0) {
					sqlite3_blob_write(blob, buff, len, offset);
					offset += len;
				} else {
					break;
				}
			}

			sqlite3_blob_close(blob);
			sqlite3_finalize(stmt);
		);
		blob = NULL;
		stmt = NULL;
		remain -= offset;
	}

	free(buff);
	fclose(fh);
	remove(temp_index);
	free(temp_index);
}

void pyfastx_load_gzip_index(char* index_file, zran_index_t* gzip_index, sqlite3* index_db) {
	int fd;
	int i, rows;

	char *temp_index;
	void *buff;

	FILE *fh;

	sqlite3_blob *blob;
	sqlite3_stmt *stmt;

	Py_ssize_t bytes = 0;
	Py_ssize_t offset;
	Py_ssize_t len;

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(index_db, "SELECT COUNT(1) FROM gzindex", -1, &stmt, NULL);
		sqlite3_step(stmt);
		rows = sqlite3_column_int(stmt, 0);
		sqlite3_finalize(stmt);
	);

	if (!rows) {
		PyErr_SetString(PyExc_RuntimeError, "no gzip index exists in index file");
		return;
	}

	temp_index = generate_random_name(index_file);
	if ((fd = mkstemp(temp_index)) < 0) {
		free(temp_index);
		PyErr_SetString(PyExc_RuntimeError, "failed to create temp file");
		return;
	}
	close(fd);

	fh = fopen(temp_index, "wb");
	buff = malloc(1048576);

	for (i = 1; i <= rows; i++) {
		offset = 0;

		PYFASTX_SQLITE_CALL(
			sqlite3_blob_open(index_db, "main", "gzindex", "content", i, 0, &blob);
			bytes = sqlite3_blob_bytes(blob);

			while (offset < bytes) {
				len = bytes - offset;

				if (len > 1048576) {
					len = 1048576;
				}

				sqlite3_blob_read(blob, buff, len, offset);
				fwrite(buff, 1, len, fh);
				offset += len;
			}
			sqlite3_blob_close(blob);
		);
		blob = NULL;
	}
	free(buff);
	fclose(fh);

	fh = fopen(temp_index, "rb");
	if (zran_import_index(gzip_index, fh, NULL) != ZRAN_IMPORT_OK) {
		PyErr_SetString(PyExc_RuntimeError, "failed to import gzip index");
	}
	fclose(fh);
	remove(temp_index);
	free(temp_index);
}

char *str_n_str(char *haystack, char *needle, Py_ssize_t len, Py_ssize_t size) {
	char *result;
	Py_ssize_t pos;

	result = strstr(haystack, needle);

	if (result) {
		pos = result - haystack + len;
		if (pos <= size) {
			return result;
		}
	}

	return NULL;
}