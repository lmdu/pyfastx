#/bin/bash
seqkit subseq --chr $1 $2 | seqkit seq -r -v -p > /dev/null 2>&1
