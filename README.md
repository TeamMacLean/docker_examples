## README

This repo contains some examples of docker files to create docker images of some commonly used bioinformatics tools - trimmomatic, bowtie2, bwa, samtools etc

## How to create docker images and run

You should have docker installed on your system. On your ubuntu system, install docker with this command:
```
sudo apt-get install -y docker
```

Now to create Trimmomatic docker image, use this command
```
sudo docker build -t dtrimmomatic ./trimmomatic
```

The option -t is for providing tag name to the docker image. You can give any tag name.

Once the image is created, use this command to see the help options for trimmomatic

```
sudo docker run -i dtrimmomatic TrimmomaticPE
sudo docker run -i dtrimmomatic TrimmomaticSE
sudo docker run -i dtrimmomatic TrimmomaticPE -version   # to get version
sudo docker run -i dtrimmomatic TrimmomaticPE -threads 2 -phred33 --trimlog trim.log --validatePairs R1.fq.gz R1.fq.gz R1paired.fq.gz R1unpaired.fq.gz R2paired.fq.gz R2unpaired.fq.gz    # sample command
```

To create image for bowtie2, bwa aligners, use this command

```
sudo docker build -t aligner ./aligners
```

To run bowtie2 and bwa commands,
```
sudo docker run -i aligner bowtie2 --help
sudo docker run -i aligner bowtie2 --version
sudo docker run -i aligner bwa mem
```

To create image for samtools, run
```
sudo docker build -t dsamtools samtools
```

Similarly, to run samtools from docker image
```
sudo docker run -t dsamtools samtools
sudo docker run -t dsamtools samtools view sample.bam
sudo docker run -t dsamtools samtools sort -o out_sorted.bam sample.bam
```

## Creating docker image for bowtie2 alignment

To create an image

```
sudo docker build -t bowtie2align bowtie2alignment
```

To run the image
```
sudo docker run -i bowtie2align 

Bash program to align short paired end reads to a reference sequence using bowtie2
----------------------------------------------------------------------------------
Usage
		-s	sample name
		-r	reference fasta file
		-1	forward reads R1
		-2	reverse reads R2
		-o	output sorted bam filename


```

Now, run bowtie2 alignment

```
sudo docker run --rm -v /home/shrestha/Downloads:/myvol  -w /myvol -i bowtie2align bowtie2alignment -s test -r Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa -1 TestDataSet_R1.fastq -2 TestDataSet_R2.fastq -o Testalign_sorted.bam
```
Options explaination:

--rm	removes the container after the execution completes
-v	mounts the directory /home/shrestha/Downloads/ to a docker volumn named /myvol
-w	go to the working directory /myvol which is the Downloads folder we have mounted.

All our input and output files are now in /myvol so, we don't have to give fullpath to the input/output files. The output files will be in the Downloads folder, no matter whereever you are when  you run the docker command to align reads.
 



