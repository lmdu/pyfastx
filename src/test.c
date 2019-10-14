kstream *ks;
kstring_t line = {0, 0, 0};

ks = ks_init(gzfp);

while (ks_getuntil(ks, '\n', &line, 0) >= 0) {
	if (line[0] == '>') {
		
	}




}

ks_destroy(ks);
free(line.s);