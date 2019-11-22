#!/bin/sh
module load R/3.2.0

for file in `ls Trimmed/*_R1trim.fastq`
do
	sampleprefix=${file#Trimmed/}
	sampleprefix=${sampleprefix%_R1trim.fastq}
	workingdir=`pwd`/
	genome="mm10"
	insertseq="MMuLV"
    echo $sampleprefix $workingdir $genome $insertseq
    
    out=${sampleprefix}.Rout
    out="Mapped/$out"
    printf "\n $out \n"

    R CMD BATCH --vanilla --slave "--args $sampleprefix $workingdir $genome $insertseq" 2-Map.R $out
done