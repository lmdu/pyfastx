#include <zlib.h>
#include "util.h"

//check file is whether exists in disk
int file_exists(char *file_name){
	if(access(file_name, F_OK) != -1){
		return 1;
	}
	return 0;
}

void upper_string(char *str){
	int i;
	for(i=0; str[i]; i++){
		str[i] = toupper(str[i]);
	}
}

PyObject *clean_seq(PyObject *self, PyObject *args){
	char *seq;
	if (!PyArg_ParseTuple(args, "s", &seq)){
		return NULL;
	}
	int i;
	int j = 0;
	for(i=0; seq[i]; i++){
		if(!isspace(seq[i])){
			seq[j++] = toupper(seq[i]);
		}
	}
	seq[j] = '\0';
	return Py_BuildValue("s", seq);
}

void truncate_seq(char *seq, int start, int end){
	int len;
	int i = 0;
	int j = 0;
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
	int start;
	int end;

	if (!PyArg_ParseTuple(args, "sii", &seq, &start, &end)){
		return NULL;
	}
	int i;
	int j = 0;
	int real_pos = 0;
	int flag;
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

//check input file is whether gzip file
int is_gzip(FILE* fd){
	int ret;
	unsigned char magic[4] = {0};
	ret = fread(magic, 1, sizeof(magic), fd);
	rewind(fd);

	if(ret != sizeof(magic)){
		return 0;
	}
	if(magic[0] != 0x1f || magic[1] != 0x8b || magic[2] != 0x08){
		return 0;
	}
	return 1;
}