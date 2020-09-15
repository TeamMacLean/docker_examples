#!/bin/bash


usage()
{
	echo "Bash program to align short paired end reads to a reference sequence using bowtie2"
	echo "----------------------------------------------------------------------------------"
	echo "Usage
		-s	sample name
		-r	reference fasta file
		-1	forward reads R1
		-2	reverse reads R2
		-o	output sorted bam filename"
}

while getopts s:r:1:2:o:h option
do
	case "$option" in
		s) sample=$OPTARG;;
		r) reference=$OPTARG;;
		1) R1=$OPTARG;;
		2) R2=$OPTARG;;
		o) output=$OPTARG;;
		h) usage; exit 0;;
	esac

done

if [ "$1" == "" ]; then
	usage
	exit 0
fi	

echo "You supplied following :"
echo "Sample name "$sample
echo Reference fasta file $reference
echo Forward reads file $R1
echo Reverse reads file $R2
echo Output bam file name $output


echo Building the reference index now

bowtie2-build -f $reference $reference
echo Building reference index completed.

echo Running bowtie2 alignment now

bowtie2  --no-discordant  --no-unal --rg-id $sample --rg SM:$sample --rg LB:NA --rg PL:Illumina  -x $reference -1 ${R1} -2 ${R2} | samtools view -hbS | samtools sort -o ${output} && samtools index ${output}
