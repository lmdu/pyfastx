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
printf "genome\tsize\tcount\tbioperl\t\tbiopython\t\tpyfaidx\t\tpyfasta\t\tpysam\t\tpyfastx1\t\tpyfastx2\t\tpyfastx3\t\tpyfastx4\t\n"

for gfile in ${gfiles[@]:1}; do
	memorys=()
	times=()
	filename=$(basename $gfile)
	filename="${filename%.*}"

	#get longest sequence name
	array=($(python3 get_fasta_info.py $gfile))

	#genome size
	gsize=${array[0]}

	#sequence counts in genome
	seqcounts=${array[1]}

	for i in $(seq 1 $repeats); do
		num=-1

		#bioperl
		measure_memory_time "perl bioperl_fasta_sequence_iterate.pl $gfile"

		#biopython
		measure_memory_time "python3 biopython_fasta_sequence_iterate.py $gfile"
		
		#pyfaidx
		measure_memory_time "python3 pyfaidx_fasta_sequence_iterate.py $gfile"

		#pyfasta
		measure_memory_time "python3 pyfasta_fasta_sequence_iterate.py $gfile"

		#pysam1 with index
		#measure_memory_time "python3 pysam_fasta_sequence_iterate_with_index.py $gfile"

		#pysam2 without index
		measure_memory_time "python3 pysam_fasta_sequence_iterate.py $gfile"

		#pyfastx1 with index
		measure_memory_time "python3 pyfastx_fasta_sequence_iterate_with_index.py $gfile"

		#pyfastx2 without index
		measure_memory_time "python3 pyfastx_fasta_sequence_iterate_without_index.py $gfile"

		#pyfastx3 with gzip index
		measure_memory_time "python3 pyfastx_fasta_sequence_iterate_with_index.py $gfile.gz"

		#pyfastx4 without gzip index
		measure_memory_time "python3 pyfastx_fasta_sequence_iterate_without_index.py $gfile.gz"
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
