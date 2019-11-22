#!/bin/sh
module load R/3.2.0
module load fastqc

mkdir FastQC
echo 'FastQC'
fastqc Reads/*_R1_001.fastq.gz --outdir=FastQC
fastqc Reads/*_R2_001.fastq.gz --outdir=FastQC

mkdir Trimmed
mkdir Temp
mkdir Mapped
mkdir Tables

for file in `ls Reads/*_R1_001.fastq.gz`
do
	sampleprefix=${file#Reads/}
	sampleprefix=${sampleprefix%_R1_001.fastq.gz}
	workingdir=`pwd`/
    echo $sampleprefix $workingdir
    
    out=${sampleprefix}.Rout
    out="Trimmed/$out"
    printf "\n $out \n"

    R CMD BATCH --vanilla --slave "--args $sampleprefix $workingdir" 1-Trim.R $out
done

mkdir Trimmed/FastQC
echo 'FastQC'
fastqc Trimmed/*_R1trim.fastq --outdir=Trimmed/FastQC
fastqc Trimmed/*_R2trim.fastq --outdir=Trimmed/FastQC