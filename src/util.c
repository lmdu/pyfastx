#include "util.h"

//check file is whether exists in disk
uint16_t file_exists(char *file_name){
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

/*
remove space using jump table
Reference:
https://github.com/lemire/despacer/blob/master/include/despacer.h
*/
const uint8_t jump_table[128] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
};

void remove_space(char *str, uint32_t len) {
	uint32_t i = 0, j = 0;
	while (i < len) {
		const char c = str[i++];
		str[j] = c;
		j += jump_table[(unsigned char)c];
	}
	str[j] = '\0';
}

void remove_space_uppercase(char *str, uint32_t len) {
	uint32_t i = 0, j = 0;
	while (i < len) {
		const char c = str[i++];
		str[j] = Py_TOUPPER(Py_CHARMASK(c));
		j += jump_table[(unsigned char)c];
	}
	str[j] = '\0';
}

/*void remove_space(char *str) {
	uint32_t i, j = 0;
	for(i = 0; str[i]; i++){
		if(!Py_ISSPACE(str[i])){
			str[j++] = str[i];
		}
	}
	str[j] = '\0';
}

void remove_space_uppercase(char *str) {
	uint32_t i, j = 0;
	for(i = 0; str[i]; i++){
		if(!Py_ISSPACE(str[i])){
			str[j++] = Py_TOUPPER(Py_CHARMASK(str[i]));
		}
	}
	str[j] = '\0';
}*/

void upper_string(char *str, uint32_t len) {
	uint32_t i;
	for (i = 0; i < len; i++) {
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
const uint8_t comp_map[128] = {
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
	char *p1 = seq;
	char *p2 = seq + strlen(seq) - 1;

	while (p1 <= p2) {
		const char c = comp_map[Py_CHARMASK(*p1)];
		*p1++ = comp_map[Py_CHARMASK(*p2)];
		*p2-- = c;
	}
}

void reverse_seq(char *seq) {
	char *p1 = seq;
	char *p2 = seq + strlen(seq) - 1;
	int c;

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

uint32_t sum_array(uint32_t arr[], int num) {
	int i, sum=0;
	for (i = 0; i < num; i++) {
		sum += arr[i];
	}

	return sum;
}

/*
char *int_to_str(int c) {
	char *str = (char *)malloc(2);
	str[0] = c;
	str[1] = '\0';
	return str;
}*/

int is_subset(char *seq1, char *seq2) {
	int i, j, m, n;
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

/*
PyObject *clean_seq(PyObject *self, PyObject *args){
	char *seq;
	if (!PyArg_ParseTuple(args, "s", &seq)){
		return NULL;
	}
	uint32_t i;
	uint32_t j = 0;
	for(i=0; seq[i]; i++){
		if(!isspace(seq[i])){
			seq[j++] = toupper(seq[i]);
		}
	}
	seq[j] = '\0';
	return Py_BuildValue("s", seq);
}

void truncate_seq(char *seq, uint32_t start, uint32_t end){
	uint32_t len;
	uint32_t i = 0;
	uint32_t j = 0;
	len = end - start + 1;

	for(i=0; i<strlen(seq); i++){
		if(isspace(seq[i])){
			continue;
		}

		seq[j++] = toupper(seq[i]);

		if(j > len){
			break;
		}
	}
	seq[j] = '\0';
}


PyObject *sub_seq(PyObject *self, PyObject *args){
	char *seq;
	uint32_t start;
	uint32_t end;

	if (!PyArg_ParseTuple(args, "sii", &seq, &start, &end)){
		return NULL;
	}
	uint32_t i;
	uint32_t j = 0;
	uint32_t real_pos = 0;
	uint16_t flag;
	for(i=0; seq[i]; i++){
		flag = isspace(seq[i]);

		if(!flag){
			real_pos++;
		}

		if(real_pos > end){
			break;
		}

		if(real_pos >= start){
			if(!flag){
				seq[j++] = toupper(seq[i]);
			}
		}
	}
	seq[j] = '\0';
	return Py_BuildValue("s", seq);
}
*/

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

// read line from file
// extracted from https://github.com/lattera/freebsd/blob/master/contrib/file/getline.c
/*ssize_t get_until_delim(char **buf, int delimiter, FILE *fp) {
	char *ptr, *eptr;

	size_t bufsiz = 100;

	if (*buf == NULL) {
		if ((*buf = malloc(bufsiz)) == NULL)
			return -1;
	}

	for (ptr = *buf, eptr = *buf + bufsiz;;) {
		int c = fgetc(fp);
		if (c == -1) {
			if (feof(fp))
				return ptr == *buf ? -1 : ptr - *buf;
			else
				return -1;
		}
		*ptr++ = c;
		if (c == delimiter) {
			*ptr = '\0';
			return ptr - *buf;
		}
		if (ptr + 2 >= eptr) {
			char *nbuf;
			size_t nbufsiz = bufsiz * 2;
			ssize_t d = ptr - *buf;
			if ((nbuf = realloc(*buf, nbufsiz)) == NULL)
				return -1;
			*buf = nbuf;
			bufsiz = nbufsiz;
			eptr = nbuf + nbufsiz;
			ptr = nbuf + d;
		}
	}
}

ssize_t get_line(char **buf, FILE *fp) {
	return get_until_delim(buf, '\n', fp);
}*/

//generate random file name according to given length
char *generate_random_name(int len) {
	int i;
	const char *chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	srand((unsigned int)time(NULL));
	char *result = (char *)malloc(sizeof(char)*(len+1));
	for (i = 0; i < len; i++) {
		result[i] = chars[rand() % 62];
	}
	result[len] = '\0';
	return result;
}


void pyfastx_build_gzip_index(zran_index_t* gzip_index, sqlite3* index_db) {
	sqlite3_stmt *stmt;
	char *temp_index;
	FILE* fd;
	uint32_t fsize;
	char *buff;
	int ret;

	ret = zran_build_index(gzip_index, 0, 0);

	if (ret != 0) {
		PyErr_SetString(PyExc_RuntimeError, "Failed to build gzip index");
		return;
	}

	//create temp gzip index file
	//temp_index = (char *)malloc(strlen(index_file) + 5);
	//strcpy(temp_index, index_file);
	//strcat(temp_index, ".tmp");
	temp_index = generate_random_name(20);
	
	fd = fopen(temp_index, "wb+");
	
	if (zran_export_index(gzip_index, fd) != ZRAN_EXPORT_OK){
		fclose(fd);
		PyErr_SetString(PyExc_RuntimeError, "Failed to export gzip index.");
		return;
	}
	
	fsize = ftell(fd);
	rewind(fd);

	buff = (char *)malloc(fsize + 1);

	if (fread(buff, fsize, 1, fd) != 1) {
		free(buff);
		return;
	}

	buff[fsize] = '\0';
	
	fclose(fd);
	remove(temp_index);
	free(temp_index);

	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(index_db, "INSERT INTO gzindex VALUES (?,?)", -1, &stmt, NULL);
		sqlite3_bind_null(stmt, 1);
		sqlite3_bind_blob(stmt, 2, buff, fsize, NULL);
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	);
	free(buff);
}

void pyfastx_load_gzip_index(zran_index_t* gzip_index, sqlite3* index_db) {
	sqlite3_stmt *stmt;
	uint32_t bytes = 0;
	FILE *fh;
	int ret;
	char *temp_index;
	
	PYFASTX_SQLITE_CALL(
		sqlite3_prepare_v2(index_db, "SELECT content FROM gzindex;", -1, &stmt, NULL);
		ret = sqlite3_step(stmt);
	);

	if (ret == SQLITE_ROW){
		PYFASTX_SQLITE_CALL(bytes = sqlite3_column_bytes(stmt, 0));
	} else {
		PYFASTX_SQLITE_CALL(sqlite3_finalize(stmt));
		PyErr_SetString(PyExc_RuntimeError, "failed to get bytes of index");
		return;
	}

	/*temp_index = (char *)malloc(strlen(index_file) + 5);
	strcpy(temp_index, index_file);
	strcat(temp_index, ".tmp");*/
	temp_index = generate_random_name(20);

	fh = fopen(temp_index, "wb");
	PYFASTX_SQLITE_CALL(
		fwrite(sqlite3_column_blob(stmt, 0), bytes, 1, fh);
		sqlite3_finalize(stmt);
	);
	fclose(fh);

	fh = fopen(temp_index, "rb");
	if (zran_import_index(gzip_index, fh) != ZRAN_IMPORT_OK) {
		fclose(fh);
		PyErr_SetString(PyExc_RuntimeError, "failed to import gzip index");
		return;
	}
	fclose(fh);
	remove(temp_index);
	free(temp_index);
}

//return large string and release memory
/*PyObject* make_large_sequence(char *seq) {
	PyObject *obj = Py_BuildValue("s", seq);
	free(seq);
	return obj;
}*/

//integer check and coversion
/*int integer_check(PyObject* num) {
	if (PyInt_Check(num) || PyLong_Check(num)) {
		return 1;
	}

	return 0;
}

int64_t integer_to_long(PyObject* num) {
	if (PyInt_Check(num)) {
		return PyInt_AsLong(num);
	} else if (PyLong_Check(num)) {
		return PyLong_AsLong(num);
	}

	PyErr_SetString(PyExc_ValueError, "the object is not an integer");
	return 0;
}*/
