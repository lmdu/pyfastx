/* The MIT License

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Last Modified: 05MAR2012 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "kseq.h"

kstream_t *ks_init(gzFile f)						
{																
	kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	
	ks->f = f;													
	ks->buf = (unsigned char*)malloc(BUF_SIZE);				
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
	if (ks_err(ks)) return -3;							
	if (ks->is_eof && ks->begin >= ks->end) return -1;	
	if (ks->begin >= ks->end) {							
		ks->begin = 0;									
		ks->end = gzread(ks->f, ks->buf, BUF_SIZE);	
		if (ks->end == 0) { ks->is_eof = 1; return -1;}	
		if (ks->end == -1) { ks->is_eof = 1; return -3;}
	}													
	return (int)ks->buf[ks->begin++];					
}

Py_ssize_t ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append)
{																	
	int gotany = 0;													
	if (dret) *dret = 0;											
	str->l = append? str->l : 0;									
	for (;;) {														
		Py_ssize_t i;														
		if (ks_err(ks)) return -3;									
		if (ks->begin >= ks->end) {									
			if (!ks->is_eof) {										
				ks->begin = 0;										
				ks->end = gzread(ks->f, ks->buf, BUF_SIZE);		
				if (ks->end == 0) { ks->is_eof = 1; break; }		
				if (ks->end == -1) { ks->is_eof = 1; return -3; }	
			} else break;											
		}															
		if (delimiter == KS_SEP_LINE) { 
			unsigned char *sep = (unsigned char*)memchr(ks->buf + ks->begin, '\n', ks->end - ks->begin);
			i = sep != NULL ? sep - ks->buf : ks->end;
		} else if (delimiter > KS_SEP_MAX) {						
			for (i = ks->begin; i < ks->end; ++i)					
				if (ks->buf[i] == delimiter) break;					
		} else if (delimiter == KS_SEP_SPACE) {						
			for (i = ks->begin; i < ks->end; ++i)					
				if (isspace(ks->buf[i])) break;						
		} else if (delimiter == KS_SEP_TAB) {						
			for (i = ks->begin; i < ks->end; ++i)					
				if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; 
		} else i = 0; /* never come to here! */						
		if (str->m - str->l < (i - ks->begin + 1)) {		
			str->m = str->l + (i - ks->begin) + 1;					
			kroundup32(str->m);										
			str->s = (char*)realloc(str->s, str->m);				
		}															
		gotany = 1;													
		memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); 
		str->l = str->l + (i - ks->begin);							
		ks->begin = i + 1;											
		if (i < ks->end) {											
			if (dret) *dret = ks->buf[i];							
			break;													
		}															
	}																
	if (!gotany && ks_eof(ks)) return -1;							
	if (str->s == 0) {												
		str->m = 1;													
		str->s = (char*)calloc(1, 1);								
	} else if (delimiter == KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; 
	str->s[str->l] = '\0';											
	return str->l;													
} 

Py_ssize_t ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) 
{ return ks_getuntil2(ks, delimiter, str, dret, 0); }

void kseq_rewind(kseq_t *ks)
{ (ks)->last_char = (ks)->f->is_eof = (ks)->f->begin = (ks)->f->end = 0; }

kseq_t *kseq_init(gzFile fd)									
{																
	kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					
	s->f = ks_init(fd);												
	return s;														
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
   -3   error reading stream
 */

Py_ssize_t kseq_read(kseq_t *seq) 
{ 
	int c;
	Py_ssize_t r; 
	kstream_t *ks = seq->f; 
	if (seq->last_char == 0) { /* then jump to the next header line */ 
		while ((c = ks_getc(ks)) >= 0 && c != '>' && c != '@'); 
		if (c < 0) return c; /* end of file or error*/ 
		seq->last_char = c; 
	} /* else: the first header char has been read in the previous call */ 
	seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ 
	if ((r=ks_getuntil(ks, 0, &seq->name, &c)) < 0) return r;  /* normal exit: EOF or error */ 
	if (c != '\n') ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */ 
	if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */ 
		seq->seq.m = 256; 
		seq->seq.s = (char*)malloc(seq->seq.m); 
	} 
	while ((c = ks_getc(ks)) >= 0 && c != '>' && c != '+' && c != '@') { 
		if (c == '\n') continue; /* skip empty lines */ 
		seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ 
		ks_getuntil2(ks, KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */ 
	} 
	if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	
	if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ 
		seq->seq.m = seq->seq.l + 2; 
		kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ 
		seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); 
	} 
	seq->seq.s[seq->seq.l] = 0;	/* null terminated string */ 
	if (c != '+') return seq->seq.l; /* FASTA */ 
	if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ 
		seq->qual.m = seq->seq.m; 
		seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); 
	} 
	while ((c = ks_getc(ks)) >= 0 && c != '\n'); /* skip the rest of '+' line */ 
	if (c == -1) return -2; /* error: no quality string */ 
	while ((c = ks_getuntil2(ks, KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l)); 
	if (c == -3) return -3; /* stream error */ 
	seq->last_char = 0;	/* we have not come to the next header line */ 
	if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ 
	return seq->seq.l; 
}

