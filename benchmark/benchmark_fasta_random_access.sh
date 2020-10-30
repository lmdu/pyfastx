#!/bin/bash

#store benchmark time and memory
tempfile=time_mem.tmp

#record memory
memorys=()

#record elapsed time
times=()

#number of programs
num=-1

#number of repeat tests
repeats=$1

#input fasta files
gfiles=$@

measure_memory_time(){
	/usr/bin/time -f "%e %M" -o $tempfile $1 > /dev/null 2>&1

	let num++

	if [ ! ${memorys[$num]} > 0 ]; then
		memorys[$num]=0
		times[$num]=0
	fi

	arr=($(cat $tempfile))

	#clear temp file
	if [ -e "$tempfile" ]; then
		rm "$tempfile"
	fi

	times[$num]=$(echo "${arr[0]}+${times[$num]}" | bc)
	memorys[$num]=$(echo "${arr[1]}+${memorys[$num]}" | bc)
}

#print header
printf "genome\tsize\tcount\tbioperl\t\tbiopython\t\tpyfaidx\t\tpyfasta\t\tpysam\t\tsamtools\t\tseqkit\t\tpyfastx\t\tpyfastx_gzip\t\n"

for gfile in ${gfiles[@]:1}; do
	memorys=()
	times=()
	filename=$(basename $gfile)
	filename="${filename%.*}"

	#get genome information
	array=($(python3 get_fasta_info.py $gfile))
	
	#genome size
	gsize=${array[0]}

	#sequence counts in genome
	seqcounts=${array[1]}

	for i in $(seq 1 $repeats); do
		num=-1

		#bioperl
		measure_memory_time "perl bioperl_fasta_random_access.pl $gfile.list $gfile"

		#biopython
		measure_memory_time "python3 biopython_fasta_random_access.py $gfile.list $gfile"
		
		#pyfaidx
		measure_memory_time "python3 pyfaidx_fasta_random_access.py $gfile.list $gfile"

		#pyfasta
		measure_memory_time "python3 pyfasta_fasta_random_access.py $gfile.list $gfile"

		#pysam
		measure_memory_time "python3 pysam_fasta_random_access.py $gfile.list $gfile"

		#samtools
		measure_memory_time "samtools faidx -r $gfile.list $gfile"

		#seqkit
		measure_memory_time "seqkit faidx --infile-list $gfile.list $gfile"

		#pyfastx
		measure_memory_time "python3 pyfastx_fasta_random_access.py $gfile.list $gfile"

		#pyfastx gzip
		measure_memory_time "python3 pyfastx_fasta_random_access.py $gfile.list $gfile.gz"
	done

	#print result
	printf "%s\t%s\t%s" $filename $gsize $seqcounts
	for((i=0;i<=$num;i++)); do
		mm=$(echo "${memorys[$i]}/$repeats" | bc)
		et=$(echo "${times[$i]}/$repeats" | bc)
		printf "\t%d\t%d" $mm $et
	done
	printf "\n"
done
