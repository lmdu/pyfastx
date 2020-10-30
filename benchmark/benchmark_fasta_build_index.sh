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

#print header
printf "genome\tsize\tcount\tbioperl\t\t\tbiopython\t\t\tpyfaidx\t\t\tpyfasta\t\t\tpysam\t\t\tsamtools\t\t\tseqkit\t\t\tpyfastx\t\t\tpyfastx_gzip\t\t\n"

for gfile in ${gfiles[@]:1}; do
	memorys=()
	times=()
	filename=$(basename $gfile)
	filename="${filename%.*}"

	for i in $(seq 1 $repeats); do
		num=-1

		#bioperl
		if [ -e "$gfile.index" ]; then
			rm "$gfile.index"
		fi
		measure_memory_time "perl bioperl_fasta_build_index.pl $gfile"
		let sizes[$num]=$(stat -c %s $gfile.index)

		#biopython
		if [ -e "$gfile.db" ]; then
			rm "$gfile.db"
		fi
		measure_memory_time "python3 biopython_fasta_build_index.py $gfile"
		let sizes[$num]=$(stat -c %s $gfile.db)

		#pyfaidx
		if [ -e "$gfile.fai" ]; then
			rm "$gfile.fai"
		fi
		measure_memory_time "python3 pyfaidx_fasta_build_index.py $gfile"
		let sizes[$num]=$(stat -c %s $gfile.fai)

		#pyfasta
		if [ -e "$gfile.gdx" ]; then
			rm "$gfile.gdx"
			rm "$gfile.flat"
		fi
		measure_memory_time "python3 pyfasta_fasta_build_index.py $gfile"
		let sizes[$num]=$(stat -c %s $gfile.gdx)

		#pysam
		if [ -e "$gfile.fai" ]; then
			rm "$gfile.fai"
		fi
		measure_memory_time "python3 pysam_fasta_build_index.py $gfile"
		let sizes[$num]=$(stat -c %s $gfile.fai)

		#samtools
		if [ -e "$gfile.fai" ]; then
			rm "$gfile.fai"
		fi
		measure_memory_time "samtools faidx $gfile"
		let sizes[$num]=$(stat -c %s $gfile.fai)

		#seqkit
		if [ -e "$gfile.seqkit.fai" ]; then
			rm "$gfile.seqkit.fai"
		fi
		measure_memory_time "seqkit faidx -f $gfile"
		let sizes[$num]=$(stat -c %s $gfile.seqkit.fai)

		#pyfastx
		if [ -e "$gfile.fxi" ]; then
			rm "$gfile.fxi"
		fi
		measure_memory_time "python3 pyfastx_fasta_build_index.py $gfile"
		let sizes[$num]=$(stat -c %s $gfile.fxi)

		#pyfastx gzip
		if [ -e "$gfile.gz.fxi" ]; then
			rm "$gfile.gz.fxi"
		fi
		measure_memory_time "python3 pyfastx_fasta_build_index.py $gfile.gz"
		let sizes[$num]=$(stat -c %s $gfile.gz.fxi)

	done

	#get genome information
	array=($(python3 get_fasta_info.py $gfile))
	
	#genome size
	gsize=${array[0]}

	#sequence counts in genome
	seqcounts=${array[1]}

	#print result
	printf "%s\t%s\t%s" $filename $gsize $seqcounts
	for((i=0;i<=$num;i++)); do
		mm=$(echo "${memorys[$i]}/$repeats" | bc)
		et=$(echo "${times[$i]}/$repeats" | bc)
		is=$(python3 -c "import math; print(math.ceil(${sizes[$i]}/1024))")
		printf "\t%d\t%d\t%d" $mm $et $is
	done
	printf "\n"
done
