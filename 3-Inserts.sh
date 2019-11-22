#!/bin/bash

# R script groups together the reads that map the the same location, counts the UMIs observed for each integration and lists them in a table and saves the output as RData files

module load R/3.2.0
mkdir Inserts

for file in `ls Tables/*_mapped_reads_table.RData`
do
	sampleprefix=${file#Tables/}
	sampleprefix=${sampleprefix%_mapped_reads_table.RData}
	workingdir=`pwd`/
    out=${sampleprefix}_inserts.Rout
    out="Inserts/$out"

    echo $sampleprefix $workingdir $out
     
    R CMD BATCH --vanilla --slave "--args $sampleprefix $workingdir" 3-Inserts.R $out
done




