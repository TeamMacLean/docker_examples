#!/bin/bash

### program to quality control short reads using trimmomatic ###
### Author: Ram Krishna Shrestha ###

usage(){


	echo """
	
	
	Program to quality control illumina short reads using trimmomatic 
	 
 	Usage
	
	 -d		output directory path. Creates if not exists
	 -o		output filename prefix
	 -1		R1 or forward reads
	 -2		R2 or reverse reads
	 -l		minium readlength to keep after trimming
       	 -q		average window quality to trim
	 -w		window size
	 -t		number of threads/processors to use
	"""
}

while getopts 1:2:l:q:w:t:d:o:h option
do
	case "$option" in 
		1)	r1=$OPTARG;;
		2)	r2=$OPTARG;;
		l)	minlen=$OPTARG;;
		q)	minqual=$OPTARG;;
		w)	windowsize=$OPTARG;;
		t)	threads=$OPTARG;;
		d)	outputdir=$OPTARG;;
		o)	outputprefix=$OPTARG;;
		h)	usage; exit  0;;
	esac
done


echo Running Trimmomatic
echo Options provided are:
echo Output directory $outputdir
echo Output filename prefix $outputprefix
echo R1	$r1
echo R2 $r2
echo minlength $minlen
echo minqual $minqual
echo window size $windowsiz
echo threads	$threads


#output_file_basename=$(basename  $r1 | sed 's/.fq$//;s/.fastq$//;s/.fq.gz$//;s/.fastq.gz//;s/_1_//;s/_R1_//;s/_1//;s/_R1//')
output_file_basename=$outputprefix

if [ ! -e $outputdir ]; then mkdir -p $outputdir; fi

cmd="TrimmomaticPE -threads $threads -phred33 -trimlog /dev/null $r1 $r2 ${outputdir}/${output_file_basename}_1_paired.fastq.gz ${outputdir}/${output_file_basename}_1_unpaired.fastq.gz ${outputdir}/${output_file_basename}_2_paired.fastq.gz ${outputdir}/${output_file_basename}_2_unpaired.fastq.gz SLIDINGWINDOW:${windowsize}:${minqual} MINLEN:${minlen}"

echo
echo Command executing:

echo $cmd

$cmd

exit 0


