#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "util.h"
#include "math.h"

#ifdef _WIN32
#include "windows.h"
#include "io.h"
/*int mkstemp(char *template) {
	if (_mktemp_s(template, strlen(template) + 1) != 0) {
		return -1;
	}

	return open(template, _O_CREAT | _O_EXCL, _S_IWRITE);
}*/
gzFile pyfastx_gzip_open(PyObject *path, const char *mode) {
	gzFile f;
	wchar_t wmode[10];
	int usize;

	if (!PyUnicode_Check(path)) {
		PyErr_Format(PyExc_TypeError,
					"str file path expected under Windows, got %R",
					Py_TYPE(path));
		return NULL;
	}

	wchar_t *wpath = PyUnicode_AsWideCharString(path, NULL);
	if (wpath == NULL) {
		return NULL;
	}

	usize = MultiByteToWideChar(CP_ACP, 0, mode, -1, wmode, Py_ARRAY_LENGTH(wmode));
	if (usize == 0) {
		PyErr_SetFromWindowsErr(0);
		PyMem_Free(wpath);
		return NULL;
	}

	Py_BEGIN_ALLOW_THREADS
	f = gzopen_w(wpath, wmode);
	Py_END_ALLOW_THREADS

	PyMem_Free(wpath);

	return f;
}

#else
#include <fcntl.h>
static uint32_t max(uint32_t a, uint32_t b) {

  if (a > b) return a;
  else       return b;
}

gzFile pyfastx_gzip_open(PyObject *path, const char *mode) {
	gzFile f;
	const char *path_str;
	PyObject *bytes;

	if (!PyUnicode_FSConverter(path, &bytes)) {
		return NULL;
	}

	path_str = PyBytes_AS_STRING(bytes);

	Py_BEGIN_ALLOW_THREADS
	f = gzopen(path_str, mode);
	Py_END_ALLOW_THREADS

	Py_DECREF(bytes);

	return f;
}
#endif

//const char ZRAN_INDEX_FILE_ID[] = {'G', 'Z', 'I', 'D', 'X'};
//const uint8_t ZRAN_INDEX_FILE_VERSION = 1;

//check file is whether exists in disk
int file_exists(PyObject *file_obj){
	FILE *file;

	if((file = _Py_fopen_obj(file_obj, "r"))){
		fclose(file);
		return 1;
	}

	PyErr_Clear();
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
int is_gzip_format(PyObject* file_obj){
	int ret;
	FILE* fd;
	unsigned char magic[4] = {0};

	fd = _Py_fopen_obj(file_obj, "rb");
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
/*char *generate_random_name(char* index_file) {
	int file_len;
	char *result;

	file_len = strlen(index_file);
	result = (char *)malloc(file_len + 8);
	sprintf(result, "%s.XXXXXX", index_file);

	return result;
}*/

/*Py_ssize_t pyfastx_gzip_index_size(zran_index_t* gzip_index) {
	zran_point_t *point;
	zran_point_t *list_end;

	uint8_t flags = 0; 

	Py_ssize_t index_size = 0;

	//ID and version size
	index_size += sizeof(ZRAN_INDEX_FILE_ID);
	index_size += 1;

	//flags size
	index_size += 1;

	//compressed size
	index_size += sizeof(gzip_index->compressed_size);

	//uncompressed size
	index_size += sizeof(gzip_index->uncompressed_size);

	//spacing size
	index_size += sizeof(gzip_index->spacing);

	//window size
	index_size += sizeof(gzip_index->window_size);

	//number of points size
	index_size += sizeof(gzip_index->npoints);

	//all points size
	point = gzip_index->list;
	list_end = gzip_index->list + gzip_index->npoints;

	while (point < list_end) {
		index_size += sizeof(point->cmp_offset);
		index_size += sizeof(point->uncmp_offset);
		index_size += sizeof(point->bits);
		flags = (point->data != NULL) ? 1 : 0;
		index_size += 1;
		point++;
	}

	//window data for every point
	point = gzip_index->list;
	list_end = gzip_index->list + gzip_index->npoints;
	while (point < list_end) {
		if (point->data == NULL) {
			point++;
			continue;
		}

		index_size += gzip_index->window_size;
		point++;
	}

	return index_size;
}*/

int pyfastx_gzip_index_write(sqlite3_stmt* stmt, const void *buff, size_t bytes) {
	int ret;

	PYFASTX_SQLITE_CALL(
		ret = sqlite3_bind_null(stmt, 1);
		if (ret != SQLITE_OK) goto fail;

		ret = sqlite3_bind_blob(stmt, 2, buff, bytes, NULL);
		if (ret != SQLITE_OK) goto fail;

		ret = sqlite3_step(stmt);
		if (ret != SQLITE_DONE) goto fail;

		ret = sqlite3_reset(stmt);
		if (ret != SQLITE_OK) goto fail;
	);

	return SQLITE_OK;

fail:
	return SQLITE_ERROR;
}

int pyfastx_gzip_index_read(sqlite3_stmt* stmt, void *buff) {
	int ret;
	size_t bytes;
	
	const void* ptr;

	PYFASTX_SQLITE_CALL(
		ret = sqlite3_step(stmt);
		if (ret != SQLITE_ROW) goto fail;
		ptr = sqlite3_column_blob(stmt, 1);
		bytes = sqlite3_column_bytes(stmt, 1);
	);

	memcpy(buff, ptr, bytes);

	return SQLITE_OK;

fail:
	return SQLITE_ERROR;
}

int pyfastx_gzip_index_export(zran_index_t* gzip_index, sqlite3* index_db) {
	int ret;
	uint8_t flags = 0;

	zran_point_t *point;
	zran_point_t *list_end;

	sqlite3_stmt *stmt;

	char *sql = "PRAGMA synchronous=OFF; BEGIN TRANSACTION;";
	PYFASTX_SQLITE_CALL(ret = sqlite3_exec(index_db, sql, NULL, NULL, NULL));
	if (ret != SQLITE_OK) goto fail;

	PYFASTX_SQLITE_CALL(
		ret = sqlite3_prepare_v2(index_db, "INSERT INTO gzindex VALUES (?,?)", -1, &stmt, NULL);
	);	
	if (ret != SQLITE_OK) goto fail;

	//write ID and version
	ret = pyfastx_gzip_index_write(stmt, ZRAN_INDEX_FILE_ID, sizeof(char)*5);
	if (ret != SQLITE_OK) goto fail;
	
	ret = pyfastx_gzip_index_write(stmt, &ZRAN_INDEX_FILE_VERSION, 1);
	if (ret != SQLITE_OK) goto fail;

	//write flags
	ret = pyfastx_gzip_index_write(stmt, &flags, 1);
	if (ret != SQLITE_OK) goto fail;

	//write compressed size
	ret = pyfastx_gzip_index_write(stmt, &gzip_index->compressed_size, sizeof(gzip_index->compressed_size));
	if (ret != SQLITE_OK) goto fail;

	//write uncompressed size
	ret = pyfastx_gzip_index_write(stmt, &gzip_index->uncompressed_size, sizeof(gzip_index->uncompressed_size));
	if (ret != SQLITE_OK) goto fail;

	//write spacing
	ret = pyfastx_gzip_index_write(stmt, &gzip_index->spacing, sizeof(gzip_index->spacing));
	if (ret != SQLITE_OK) goto fail;

	//write window size
	ret = pyfastx_gzip_index_write(stmt, &gzip_index->window_size, sizeof(gzip_index->window_size));
	if (ret != SQLITE_OK) goto fail;

	//write number of points
	ret = pyfastx_gzip_index_write(stmt, &gzip_index->npoints, sizeof(gzip_index->npoints));
	if (ret != SQLITE_OK) goto fail;

	//all points offset
	point = gzip_index->list;
	list_end = gzip_index->list + gzip_index->npoints;
	while (point < list_end) {
		//write compressed offset
		ret = pyfastx_gzip_index_write(stmt, &point->cmp_offset, sizeof(point->cmp_offset));
		if (ret != SQLITE_OK) goto fail;

		//write uncompressed offset
		ret = pyfastx_gzip_index_write(stmt, &point->uncmp_offset, sizeof(point->uncmp_offset));
		if (ret != SQLITE_OK) goto fail;

		//write bit offset
		ret = pyfastx_gzip_index_write(stmt, &point->bits, sizeof(point->bits));
		if (ret != SQLITE_OK) goto fail;

		//write data flag
		flags = (point->data != NULL) ? 1 : 0;
		ret = pyfastx_gzip_index_write(stmt, &flags, 1);
		if (ret != SQLITE_OK) goto fail;

		++point;
	}

	//write window data
	point = gzip_index->list;
	list_end = gzip_index->list + gzip_index->npoints;
	while (point < list_end) {
		if (point->data == NULL) {
			++point;
			continue;
		}

		//write checkpoint data
		ret = pyfastx_gzip_index_write(stmt, point->data, gzip_index->window_size);
		if (ret != SQLITE_OK) goto fail;

		++point;
	}

	PYFASTX_SQLITE_CALL(ret = sqlite3_finalize(stmt));
	if (ret != SQLITE_OK) goto fail;

	PYFASTX_SQLITE_CALL(sqlite3_exec(index_db, "COMMIT;", NULL, NULL, NULL));

	return ZRAN_EXPORT_OK;

fail:
	return ZRAN_EXPORT_WRITE_ERROR;
}

int pyfastx_gzip_index_import(zran_index_t* gzip_index, sqlite3* index_db) {
	int ret;

	uint64_t i;
	zran_point_t *point;
	zran_point_t *list_end;

	uint8_t *dataflags = NULL;

	char file_id[sizeof(char)*5];
	uint8_t version;
	uint8_t flags;

	uint64_t compressed_size;
	uint64_t uncompressed_size;
	uint32_t spacing;
	uint32_t window_size;
	uint32_t npoints;
	zran_point_t *new_list = NULL;

	sqlite3_stmt *stmt;

	PYFASTX_SQLITE_CALL(
		ret = sqlite3_prepare_v2(index_db, "SELECT * FROM gzindex", -1, &stmt, NULL);
	);
	if (ret != SQLITE_OK) goto fail;

	gzip_index->flags |= ZRAN_SKIP_CRC_CHECK;
	
	//read and verify ID
	ret = pyfastx_gzip_index_read(stmt, file_id);
	if (ret != SQLITE_OK) goto read_error;

	if (memcmp(file_id, ZRAN_INDEX_FILE_ID, sizeof(file_id))) goto unknown_format;

	//read format version and check
	ret = pyfastx_gzip_index_read(stmt, &version);
	if (ret != SQLITE_OK) goto read_error;

	if (version > ZRAN_INDEX_FILE_VERSION) goto unsupported_version;

	//read flags
	ret = pyfastx_gzip_index_read(stmt, &flags);
	if (ret != SQLITE_OK) goto read_error;

	//read compressed size and check
	ret = pyfastx_gzip_index_read(stmt, &compressed_size);
	if (ret != SQLITE_OK) goto read_error;

	if (compressed_size != gzip_index->compressed_size) goto inconsistent;

	//read uncompressed size and check
	ret = pyfastx_gzip_index_read(stmt, &uncompressed_size);
	if (ret != SQLITE_OK) goto read_error;

	if (uncompressed_size != 0 && gzip_index->uncompressed_size != 0 && gzip_index->uncompressed_size != uncompressed_size) goto inconsistent;

	//read spacing
	ret = pyfastx_gzip_index_read(stmt, &spacing);
	if (ret != SQLITE_OK) goto read_error;

	//read window size
	ret = pyfastx_gzip_index_read(stmt, &window_size);
	if (ret != SQLITE_OK) goto read_error;

	//check spacing and window size
	if (window_size < 32768) goto fail;
	if (spacing < window_size) goto fail;

	//read no. of points
	ret = pyfastx_gzip_index_read(stmt, &npoints);
	if (ret != SQLITE_OK) goto read_error;

	new_list = calloc(1, sizeof(zran_point_t) * max(npoints, 8));
	if (new_list == NULL) goto memory_error;

	dataflags = calloc(npoints, 1);
	if (dataflags == NULL) goto memory_error;

	for (i = 0, point = new_list; i < npoints; ++i, ++point) {
		ret = pyfastx_gzip_index_read(stmt, &point->cmp_offset);
		if (ret != SQLITE_OK) goto read_error;

		ret = pyfastx_gzip_index_read(stmt, &point->uncmp_offset);
		if (ret != SQLITE_OK) goto read_error;

		ret = pyfastx_gzip_index_read(stmt, &point->bits);
		if (ret != SQLITE_OK) goto read_error;

		if (version >= 1) {
			ret = pyfastx_gzip_index_read(stmt, &flags);
			if (ret != SQLITE_OK) goto read_error;
		} else {
			flags = (point == new_list) ? 0 : 1;
		}

		dataflags[i] = flags;
	}

	for (i = 0, point = new_list; i < npoints; ++i, ++point) {
		if (dataflags[i] == 0) {
			continue;
		}

		point->data = calloc(1, window_size);
		if (point->data == NULL) goto memory_error;

		ret = pyfastx_gzip_index_read(stmt, point->data);
		if (ret != SQLITE_OK) goto read_error;
	}

	PYFASTX_SQLITE_CALL(ret = sqlite3_finalize(stmt));
	if (ret != SQLITE_OK) goto fail;

	if (gzip_index->uncompressed_size == 0 && uncompressed_size != 0) {
		gzip_index->uncompressed_size = uncompressed_size;
	}

	if (gzip_index->spacing != spacing) {
		gzip_index->spacing = spacing;
	}

	if (gzip_index->window_size != window_size) {
		gzip_index->window_size = window_size;
	}

	point = gzip_index->list + 1;
	list_end = gzip_index->list + gzip_index->npoints;

	while (point < list_end) {
		free(point->data);
		++point;
	}

	free(gzip_index->list);
	gzip_index->list = new_list;
	gzip_index->npoints = npoints;
	gzip_index->size = max(npoints, 8);
	free(dataflags);

	return ZRAN_IMPORT_OK;

fail:
	ret = ZRAN_IMPORT_FAIL;
	goto cleanup;

read_error:
	ret = ZRAN_IMPORT_READ_ERROR;
	goto cleanup;

inconsistent:
	ret = ZRAN_IMPORT_INCONSISTENT;
	goto cleanup;

memory_error:
	ret = ZRAN_IMPORT_MEMORY_ERROR;
	goto cleanup;

unknown_format:
	ret = ZRAN_IMPORT_UNKNOWN_FORMAT;
	goto cleanup;

unsupported_version:
	ret = ZRAN_IMPORT_UNSUPPORTED_VERSION;
	goto cleanup;

cleanup:
	if (new_list != NULL) {
		point    = new_list + 1;
		list_end = new_list + npoints;

		while (point < list_end && point->data != NULL) {
			free(point->data);
			++point;
		}

		free(new_list);
	}

	if (dataflags != NULL) {
		free(dataflags);
	}

	return ret;
}

void pyfastx_build_gzip_index(zran_index_t* gzip_index, sqlite3* index_db) {
	int ret;

	ret = zran_build_index(gzip_index, 0, 0);
	if (ret != 0) {
		PyErr_Format(PyExc_RuntimeError, "failed to build gzip index return %d", ret);
		return;
	}

	ret = pyfastx_gzip_index_export(gzip_index, index_db);
	if (ret != ZRAN_EXPORT_OK) {
		PyErr_Format(PyExc_RuntimeError, "failed to save gzip index return %d", ret);
		return;
	}
}

void pyfastx_load_gzip_index(zran_index_t* gzip_index, sqlite3* index_db) {
	int ret;
	int rows;

	sqlite3_stmt *stmt;

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

	ret = pyfastx_gzip_index_import(gzip_index, index_db);
	if (ret != ZRAN_IMPORT_OK) {
		PyErr_Format(PyExc_RuntimeError, "failed to import gzip index return %d", ret);
		return;
	}
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