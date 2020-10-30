#!/bin/bash

#store benchmark time and memory
tempfile=time_mem.tmp

#record memory
memorys=()

#record elapsed time
times=()

#record index file size
sizes=()

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

printf "file\tbase\tcount\tfsize\tgsize\tbioperl\t\tbiopython1\t\tbiopython2\t\tpysam\t\tpyfastx1\t\tpyfastx2\t\tpyfastx3\t\tpyfastx4\t\n"

for gfile in ${gfiles[@]:1}; do
	memorys=()
	times=()
	filename=$(basename $gfile)
	filename="${filename%.*}"

	for i in $(seq 1 $repeats); do
		num=-1

		#bioperl
		measure_memory_time "perl bioperl_fastq_sequence_iterate.pl $gfile"

		#biopython1
		measure_memory_time "python3 biopython_fastq_sequence_iterate_with_seqio.py $gfile"

		#biopython2
		measure_memory_time "python3 biopython_fastq_sequence_iterate_with_qualityio.py $gfile"

		#pysam
		measure_memory_time "python3 pysam_fastq_sequence_iterate.py $gfile"

		#pyfastx1 with index
		measure_memory_time "python3 pyfastx_fastq_sequence_iterate_with_index.py $gfile"

		#pyfastx2 without index
		measure_memory_time "python3 pyfastx_fastq_sequence_iterate_without_index.py $gfile"

		#pyfastx3 gzip with index
		measure_memory_time "python3 pyfastx_fastq_sequence_iterate_with_index.py $gfile.gz"

		#pyfastx4 gzip without index
		measure_memory_time "python3 pyfastx_fastq_sequence_iterate_without_index.py $gfile.gz"

	done

	#get fastq information
	array=($(python3 get_fastq_info.py $gfile))

	#fastq nucleotides size
	gsize=${array[0]}

	#reads counts in fastq
	seqcounts=${array[1]}

	#fastq file size
	fsize=$(stat -c %s $gfile)

	#fastq gzip file size
	gzsize=$(stat -c %s $gfile.gz)

	#print result
	printf "%s\t%s\t%s\t%s\t%s" $filename $gsize $seqcounts $fsize $gzsize
	for((i=0;i<=$num;i++)); do
		mm=$(echo "${memorys[$i]}/$repeats" | bc)
		et=$(echo "${times[$i]}/$repeats" | bc)
		#bi=$(echo "scale=2;${sizes[$i]}/1024" | bc)
		printf "\t%d\t%d" $mm $et
	done
	printf "\n"
done
