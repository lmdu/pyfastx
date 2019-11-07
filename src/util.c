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

void remove_space(char *str) {
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
			str[j++] = Py_TOUPPER((unsigned char) str[i]);
		}
	}
	str[j] = '\0';
}

void upper_string(char *str) {
	uint32_t i;
	for(i = 0; str[i]; i++) {
		if(str[i] >= 'a' && str[i] <= 'z'){
			str[i] = str[i] - 32;
		}
	}
}

void reverse_seq(char *seq) {
	char *p1 = seq;
	char *p2 = seq + strlen(seq) - 1;
	int tmp;

	while (p1 < p2) {
		tmp = *p1;
		*p1++ = *p2;
		*p2-- = tmp;
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
*/
void generate_complement_map(int arr[]) {
	arr[65] = arr[97] = 84;
	arr[71] = arr[103] = 67;
	arr[67] = arr[99] = 71;
	arr[84] = arr[116] = 65;
	arr[78] = arr[110] = 78;
	arr[77] = arr[109] = 75;
	arr[82] = arr[114] = 89;
	arr[87] = arr[119] = 87;
	arr[83] = arr[115] = 83;
	arr[89] = arr[121] = 82;
	arr[75] = arr[107] = 77;
	arr[86] = arr[118] = 66;
	arr[72] = arr[104] = 68;
	arr[68] = arr[100] = 72;
	arr[66] = arr[98] = 86;
	arr[85] = arr[117] = 65;
}

void reverse_complement_seq(char *seq) {
	char *p1 = seq;
	char *p2 = seq + strlen(seq) - 1;
	int tmp;

	//generate complement mapping
	int mapping[125];
	generate_complement_map(mapping);

	while (p1 < p2) {
		tmp = mapping[Py_CHARMASK(*p1)];
		*p1++ = mapping[Py_CHARMASK(*p2)];
		*p2-- = tmp;
	}

	if (p1 == p2) {
		*p1 = mapping[Py_CHARMASK(*p1)];
	}
}

void complement_seq(char *seq) {
	//generate complement mapping
	int mapping[125];
	int c;
	generate_complement_map(mapping);

	while ((c=*seq)) {
		*seq++ = mapping[Py_CHARMASK(c)];
	}
}

uint32_t sum_array(uint32_t arr[], int num) {
	int i, sum=0;
	for (i = 0; i < num; i++) {
		sum += arr[i];
	}

	return sum;
}

char *int_to_str(int c) {
	char *str = (char *)malloc(2);
	str[0] = c;
	str[1] = '\0';
	return str;
}

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
uint16_t is_gzip_format(char* file_name){
	uint16_t ret;
	FILE* fd;
	unsigned char magic[4] = {0};

	fd = fopen(file_name, "rb");
	ret = fread(magic, 1, sizeof(magic), fd);
	fclose(fd);

	if(ret != sizeof(magic)){
		return 0;
	}
	
	if(magic[0] != 0x1f || magic[1] != 0x8b || magic[2] != 0x08){
		return 0;
	}
	
	return 1;
}

/* read line from zran indexed gzip file
int64_t zran_readline(zran_index_t *index, char *linebuf, uint32_t bufsize) {
	int64_t startpos;
	int64_t ret;
	uint16_t haveline;
	char* readbuf = (char *)malloc(bufsize);
	char* lineidx;
	uint32_t idxpos;

	//record start position
	startpos = zran_tell(index);

	while (1) {
		ret = zran_read(index, readbuf, bufsize);

		if (ret == ZRAN_READ_EOF) {
			return 0
		}

		lineidx = strchr(readbuf, '\n');
		haveline = (lineidx != NULL) ? 1 : 0;

		if (haveline) {
			idxpos = (int)(lineidx - readbuf);
			linebuf = (char *)realloc(linebuf, strlen(linebuf) + idxpos);
			memcpy(linebuf+strlen(linebuf), readbuf, idxpos);
			linebuf[strlen[linebuf]] = '\0';
		} else {
			linebuf = (char *)realloc(linebuf, strlen(linebuf) + bufsize);
			memcpy(linebuf+strlen(linebuf), readbuf, bufsize);
			linebuf[strlen(linebuf)] = '\0';
		}

		if (haveline) {
			break;
		}
	}
	zran_seek(index, startpos+strlen(linebuf), SEEK_SET, NULL);
	free(readbuf);
	return strlen(linebuf);
}*/

void pyfastx_build_gzip_index(zran_index_t* gzip_index, sqlite3* index_db, char* index_file) {
	sqlite3_stmt *stmt;

	zran_build_index(gzip_index, 0, 0);

	//create temp gzip index file
	char *temp_index = (char *)malloc(strlen(index_file) + 5);
	strcpy(temp_index, index_file);
	strcat(temp_index, ".tmp");
	
	FILE* fd = fopen(temp_index, "wb+");
	
	if(zran_export_index(gzip_index, fd) != ZRAN_EXPORT_OK){
		PyErr_SetString(PyExc_RuntimeError, "Failed to export gzip index.");
	}
	
	uint32_t fsize = ftell(fd);
	rewind(fd);

	char *buff = (char *)malloc(fsize + 1);

	if(fread(buff, fsize, 1, fd) != 1){
		return;
	}

	buff[fsize] = '\0';
	
	fclose(fd);
	remove(temp_index);

	sqlite3_prepare_v2(index_db, "INSERT INTO gzindex VALUES (?,?)", -1, &stmt, NULL);
	sqlite3_bind_null(stmt, 1);
	sqlite3_bind_blob(stmt, 2, buff, fsize, NULL);
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
	free(buff);
}

void pyfastx_load_gzip_index(zran_index_t* gzip_index, sqlite3* index_db, char* index_file) {
	sqlite3_stmt *stmt;
	uint32_t bytes = 0;
	FILE *fh;
	
	char *temp_index = (char *)malloc(strlen(index_file) + 5);
	strcpy(temp_index, index_file);
	strcat(temp_index, ".tmp");
	
	fh = fopen(temp_index, "wb");

	sqlite3_prepare_v2(index_db, "SELECT content FROM gzindex;", -1, &stmt, NULL);
	if (sqlite3_step(stmt) == SQLITE_ROW){
		bytes = sqlite3_column_bytes(stmt, 0);
	}
	
	fwrite(sqlite3_column_blob(stmt, 0), bytes, 1, fh);
	fclose(fh);

	fh = fopen(temp_index, "rb");

	if (zran_import_index(gzip_index, fh) != ZRAN_IMPORT_OK){
		PyErr_SetString(PyExc_RuntimeError, "Failed to import gzip index.");
	}
	fclose(fh);
	remove(temp_index);
}
