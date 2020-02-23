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
	/usr/bin/time -f "%e %M" -o $tempfile $1 > /dev/null #2>&1

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

	#build index
	#if [ ! -e "$gfile.fai" ]; then
	#	samtools faidx $gfile
	#fi
	#seqkit faidx -f $gfile 2> /dev/null
	#python3 pyfastx_build_index.py $gfile
	#python3 pyfastx_build_index.py $gfile.gz
	#python3 pyfasta_build_index.py $gfile

	#get longest sequence name
	array=($(python3 get_genome_info.py $gfile))
	
	#genome size
	gsize=${array[0]}

	#sequence counts in genome
	seqcounts=${array[1]}

	#longest sequence name
	longname=${array[2]}

	#longest sequence length
	longlen=${array[3]}

	for i in $(seq 1 $repeats); do
		num=-1
		
		#pyfaidx
		measure_memory_time "python3 pyfaidx_reverse_complement.py $longname $gfile"

		#pyfasta
		measure_memory_time "python3 pyfasta_reverse_complement.py $longname $gfile"

		#samtools
		measure_memory_time "samtools faidx -i $gfile $longname"

		#seqkit
		measure_memory_time "sh seqkit_reverse_complement.sh $longname $gfile"

		#pyfastx
		measure_memory_time "python3 pyfastx_reverse_complement.py $longname $gfile"

		#pyfastx gzip
		measure_memory_time "python3 pyfastx_reverse_complement.py $longname $gfile.gz"
	done

	# clear index files
	#for idxfile in "$gfile.fxi" "$gfile.gz.fxi" "$gfile.fai" "$gfile.gdx" "$gfile.flat" "$gfile.seqkit.fai"; do
	#	if [ -e "$idxfile" ]; then
	#		rm "$idxfile"
	#	fi
	#done

	#print result
	printf "%s\t%s\t%s\t%s\t%s\t" $filename $gsize $seqcounts $longname $longlen
	for((i=0;i<=$num;i++)); do
		mm=$(echo "scale=2;${memorys[$i]}/$repeats" | bc)
		et=$(echo "scale=2;${times[$i]}/$repeats" | bc)
		printf "%.2f\t%.2f\t" $mm $et
	done
	printf "\n"
done
