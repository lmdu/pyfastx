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

for gfile in ${gfiles[@]:1}; do
	memorys=()
	times=()
	filename=$(basename $gfile)
	filename="${filename%.*}"

	for i in $(seq 1 $repeats); do
		num=-1
		
		#pyfaidx
		if [ -e "$gfile.fai" ]; then
			rm "$gfile.fai"
		fi
		measure_memory_time "python3 pyfaidx_build_index.py $gfile"

		#pyfasta
		if [ -e "$gfile.gdx" ]; then
			rm "$gfile.gdx"
			rm "$gfile.flat"
		fi
		measure_memory_time "python3 pyfasta_build_index.py $gfile"

		#pysam
		if [ -e "$gfile.fai" ]; then
			rm "$gfile.fai"
		fi
		measure_memory_time "python3 pysam_build_index.py $gfile"

		#samtools
		if [ -e "$gfile.fai" ]; then
			rm "$gfile.fai"
		fi
		measure_memory_time "samtools faidx $gfile"

		#seqkit
		if [ -e "$gfile.seqkit.fai" ]; then
			rm "$gfile.seqkit.fai"
		fi
		measure_memory_time "seqkit faidx -f $gfile"

		#pyfastx
		if [ -e "$gfile.fxi" ]; then
			rm "$gfile.fxi"
		fi
		measure_memory_time "python3 pyfastx_build_index.py $gfile"

		#pyfastx gzip
		if [ -e "$gfile.gz.fxi" ]; then
			rm "$gfile.gz.fxi"
		fi
		measure_memory_time "python3 pyfastx_build_index.py $gfile.gz"

	done

	#get genome information
	array=($(python3 get_genome_info.py $gfile))
	
	#genome size
	gsize=${array[0]}

	#sequence counts in genome
	seqcounts=${array[1]}

	#longest sequence name
	#longname=${array[2]}

	#longest sequence length
	#longlen=${array[3]}

	#print result
	printf "%s\t%s\t%s\t" $filename $gsize $seqcounts
	for((i=0;i<=num;i++)); do
		mm=$(echo "scale=2;${memorys[$i]}/$repeats" | bc)
		et=$(echo "scale=2;${times[$i]}/$repeats" | bc)
		printf "%.2f\t%.2f\t" $mm $et
	done
	printf "\n"
done
