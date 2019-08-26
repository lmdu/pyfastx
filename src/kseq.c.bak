#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "kseq.h"

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_MAX   1
#define BUFF_SIZE    40960

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)


kstream_t *ks_init(gzFile f)
{																
	kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	
	ks->f = f;													
	ks->buf = (char*)malloc(BUFF_SIZE);							
	return ks;													
}																
void ks_destroy(kstream_t *ks)					
{																
	if (ks) {													
		free(ks->buf);											
		free(ks);												
	}															
}

int ks_getc(kstream_t *ks)				
{														
	if (ks->is_eof && ks->begin >= ks->end) return -1;	
	if (ks->begin >= ks->end) {							
		ks->begin = 0;									
		ks->end = gzread(ks->f, ks->buf, BUFF_SIZE);	
		if (ks->end < BUFF_SIZE) ks->is_eof = 1;		
		if (ks->end == 0) return -1;					
	}													
	return (int)ks->buf[ks->begin++];					
}

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) 
{																	
	if (dret) *dret = 0;											
	str->l = 0;														
	if (ks->begin >= ks->end && ks->is_eof) return -1;				
	for (;;) {														
		int i;														
		if (ks->begin >= ks->end) {									
			if (!ks->is_eof) {										
				ks->begin = 0;										
				ks->end = gzread(ks->f, ks->buf, BUFF_SIZE);		
				if (ks->end < BUFF_SIZE) ks->is_eof = 1;			
				if (ks->end == 0) break;							
			} else break;											
		}															
		if (delimiter > KS_SEP_MAX) {								
			for (i = ks->begin; i < ks->end; ++i)					
				if (ks->buf[i] == delimiter) break;					
		} else if (delimiter == KS_SEP_SPACE) {						
			for (i = ks->begin; i < ks->end; ++i)					
				if (isspace(ks->buf[i])) break;						
		} else if (delimiter == KS_SEP_TAB) {						
			for (i = ks->begin; i < ks->end; ++i)					
				if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; 
		} else i = 0; /* never come to here! */						
		if (str->m - str->l < i - ks->begin + 1) {					
			str->m = str->l + (i - ks->begin) + 1;					
			kroundup32(str->m);										
			str->s = (char*)realloc(str->s, str->m);				
		}															
		memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); 
		str->l = str->l + (i - ks->begin);							
		ks->begin = i + 1;											
		if (i < ks->end) {											
			if (dret) *dret = ks->buf[i];							
			break;													
		}															
	}																
	if (str->l == 0) {												
		str->m = 1;													
		str->s = (char*)calloc(1, 1);								
	}																
	str->s[str->l] = '\0';											
	return str->l;													
}


kseq_t *kseq_init(gzFile fd)							
{																	
	kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					
	s->f = ks_init(fd);												
	return s;														
}																	
void kseq_rewind(kseq_t *ks)							
{
	ks->last_char = 0;												
	ks->f->is_eof = ks->f->begin = ks->f->end = 0;					
}
void kseq_destroy(kseq_t *ks)							
{																	
	if (!ks) return;												
	free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s); 
	ks_destroy(ks->f);												
	free(ks);														
}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
int kseq_read(kseq_t *seq)									
{																	
	int c;															
	kstream_t *ks = seq->f;											
	if (seq->last_char == 0) { /* then jump to the next header line */ 
		while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');	
		if (c == -1) return -1; /* end of file */					
		seq->last_char = c;											
	} /* the first header char has been read */						
	seq->comment.l = seq->seq.l = seq->qual.l = 0;					
	if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1;			
	if (c != '\n') ks_getuntil(ks, '\n', &seq->comment, 0);			
	while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') { 
		if (isgraph(c)) { /* printable non-space character */		
			if (seq->seq.l + 1 >= seq->seq.m) { /* double the memory */ 
				seq->seq.m = seq->seq.l + 2;						
				kroundup32(seq->seq.m); /* rounded to next closest 2^k */ 
				seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); 
			}														
			seq->seq.s[seq->seq.l++] = (char)c;						
		}															
	}																
	if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	
	seq->seq.s[seq->seq.l] = 0;	/* null terminated string */		
	if (c != '+') return seq->seq.l; /* FASTA */					
	if (seq->qual.m < seq->seq.m) {	/* allocate enough memory */	
		seq->qual.m = seq->seq.m;									
		seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m);		
	}																
	while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */ 
	if (c == -1) return -2; /* we should not stop here */			
	while ((c = ks_getc(ks)) != -1 && seq->qual.l < seq->seq.l)		
		if (c >= 33 && c <= 127) seq->qual.s[seq->qual.l++] = (unsigned char)c;	
	seq->qual.s[seq->qual.l] = 0; /* null terminated string */		
	seq->last_char = 0;	/* we have not come to the next header line */ 
	if (seq->seq.l != seq->qual.l) return -2; /* qual string is shorter than seq string */ 
	return seq->seq.l;												
}
