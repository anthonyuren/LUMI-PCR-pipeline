#######################################################################################################################
# LUMI-PCR Script for custom trimming of reads
#######################################################################################################################

# Libraries used
library(ShortRead)
library(Biostrings)
library(gtools)
sessionInfo()

#######################################################################################################################
# Input filenames/arguments
#######################################################################################################################

args <- commandArgs(trailingOnly = TRUE)
sample_prefix <- args[1] 
work_dir <- args[2]      
# e.g. sample_prefix <- "M-1_S95"; work_dir <- '/Volumes/ILLUMINA/M-1-24b/'
# e.g. sample_prefix <- "986_ATACGAAGCC-TAACAGAGAG"; work_dir <- '~/Desktop/Replicates8'
setwd(work_dir)

file1 <- paste("Reads/", sample_prefix, "_R1_001.fastq.gz", sep="") # the read 1 read e.g. Trimgalore/M-1_S95_R1_001_trimmed.fastq
file2 <- paste("Reads/", sample_prefix, "_R2_001.fastq.gz", sep="") # the read 2 read e.g. Trimgalore/M-1_S95_R2_001_trimmed.fastq
file3 <- paste("Reads/", sample_prefix, "_R3_001.fastq.gz", sep="")      # the UMI read e.g. Reads/M-1_S95_R3_001.fastq.gz

print(work_dir)

print("Read files.")
print(file1)
print(file2)
print(file3)

#######################################################################################################################
# load table of samples and adapters and the integration (virus/transposon) sequences
#######################################################################################################################

sample_table <- read.delim("sample_table.txt", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
print(colnames(sample_table))
sample_seqs <- sample_table[which(sample_table$sample_id == sample_prefix),]
print(paste("sample_id =",sample_prefix))
print(as.list(sample_seqs))

#######################################################################################################################
#######################################################################################################################
# input the read files into a single table reads_table
#######################################################################################################################

read1 <- readFastq(paste(file1))
read2 <- readFastq(paste(file2))
read3 <- readFastq(paste(file3))

# for very large read numbers from multiple flow cells can optionally process just a subset at a time
# e.g. subset <- 1:100000
# e.g. subset <- grep("HISEQ:313:HB12PADXX",as.character(id(read1)))
# read1 <- read1[subset]
# read2 <- read2[subset]
# read3 <- read3[subset]

read1_ID_split <- unlist(strsplit(as.character(id(read1)) , " "))
i <- seq(1, length(read1_ID_split), 2)
IDshort <- read1_ID_split[i]
read1_table <- as.data.frame(cbind(IDshort, as.character(id(read1)), as.character(sread(read1)), as.character(PhredQuality(quality(read1))) ), stringsAsFactors = FALSE)
names(read1_table) <- c("Id","R1_ID_long","R1_seq", "R1_qual")


read2_ID_split <- unlist(strsplit(as.character(id(read2)) , " "))
i <- seq(1, length(read2_ID_split), 2)
IDshort <- read2_ID_split[i]
read2_table <- as.data.frame(cbind(IDshort, as.character(id(read2)), as.character(sread(read2)), as.character(PhredQuality(quality(read2))) ), stringsAsFactors = FALSE)
names(read2_table) <- c("Id","R2_ID_long","R2_seq","R2_qual")


read3_ID_split <- unlist(strsplit(as.character(id(read3)) , " "))
i <- seq(1, length(read3_ID_split), 2)
IDshort <- read3_ID_split[i]
read3_table <- as.data.frame(cbind(IDshort, as.character(id(read3)), as.character(sread(read3)), as.character(PhredQuality(quality(read3))) ), stringsAsFactors = FALSE)
names(read3_table) <- c("Id","R3_ID_long","R3_seq","R3_qual")


# create a single table for all reads
read1_2_table <- merge(read1_table, read2_table, by.x=1, by.y=1, all.x=TRUE)
reads_table <- merge(read1_2_table, read3_table, by.x=1, by.y=1, all.x=TRUE)

rm(read1_ID_split, read2_ID_split, read3_ID_split, IDshort, read1_2_table, read1_table, read2_table, read3_table)

save(reads_table, file = paste("Tables/",sample_prefix,"_reads_table_full.RData",sep=""))

print(paste(nrow(reads_table),"reads in unfiltered table"))
#######################################################################################################################

#######################################################################################################################
# optional step: remove reads with an average read quality below a certain quality threshold
#######################################################################################################################

average_r1_qual <- rep(0,nrow(reads_table))
average_r2_qual <- rep(0,nrow(reads_table))
average_min_r1_r2 <- rep(0,nrow(reads_table))

for(i in 1:nrow(reads_table)) {
  average_r1_qual[i] <- mean(c(asc(reads_table$R1_qual[i]))-33)
  average_r2_qual[i] <- mean(c(asc(reads_table$R2_qual[i]))-33)
  average_min_r1_r2[i] <- min(average_r1_qual[i], average_r2_qual[i])
}

reads_table <- reads_table[(which(average_min_r1_r2 > 25)),] # threshold value
print(paste(nrow(reads_table),"reads passing quality cutoff"))
#######################################################################################################################
# find insert bases from the insert-genome junction at the start of read 2
# eliminate reads that don't have these bases and save table
#######################################################################################################################
reads_table <- cbind(reads_table,NA)
colnames(reads_table)[11] <- "R2_junction_bases"

junction_gap_len <- nchar(sample_seqs$seq_insert_bases) # how many bases from the beginning of read 2 is the junction expected

read2_first_bases <- substr(reads_table$R2_seq, start = 1, stop = junction_gap_len + 10) # record the junction + 10 genome bases (can export to view/sort/QC)
reads_table$R2_junction_bases <- read2_first_bases

reads_table <- cbind(reads_table,1)
colnames(reads_table)[12] <- "R2trim_5p"


read2_first_bases <- substr(reads_table$R2_seq, start = 1, stop = junction_gap_len)   # hist(adist(read2_start_string, seq_insert_bases))
# allow mismatches from consensus sequence since sometimes read 2 begins very poorly
read2_start_trim <- which(adist(read2_first_bases, sample_seqs$seq_insert_bases) < ceiling(sqrt(junction_gap_len)*junction_gap_len/9)   ) 
reads_table$R2trim_5p[read2_start_trim] <- junction_gap_len + 1                       #start trim at the start of the read

# limit the reads_table to only those that have expected bases of insert sequence at 5' end
reads_table  <- reads_table[which(reads_table$R2trim_5p == junction_gap_len + 1),]

save(reads_table, file = paste("Tables/",sample_prefix,"_reads_table.RData",sep=""))  # these are the only reads that will be analysed further
print(paste(nrow(reads_table),"reads passing quality cutoff with expected read 2 5' bases"))
#######################################################################################################################
# find insert sequences at end of read1
#######################################################################################################################

R1_insert_trim <- array(NA,c(nrow(reads_table),19))
R1_insert_trim <- as.data.frame(R1_insert_trim, stringsAsFactors = FALSE)
colnames(R1_insert_trim) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen","trim_coord","primer_seq","align1", "align2", "align3")

primer_seq <- sample_seqs$primer_seq_i5_r2_end # can be found at the 3' end of read 1 in reverse complement orientation
primer_seq <- paste(primer_seq, sample_seqs$pcr_insert_bases, sep="") # include bases prior to the insert junction that are not present in the secondary PCR primer
print(primer_seq)
primer_seq <- DNAStringSet(primer_seq)
writeXStringSet(primer_seq, "Temp/primer.fasta", format="fasta")

# quick loop to write a fasta file for every read1
for(read_num in 1:nrow(reads_table)) {
  fileConn <- file(paste("Temp/r1_",read_num,".fasta",sep=""))
  line1 <- paste(">",reads_table$Id[read_num],sep="")
  read_seq <- reads_table$R1_seq[read_num]
  line2 <- substr(read_seq,1,80)
  line3 <- substr(read_seq,81,160)
  line4 <- substr(read_seq,161,240)
  line5 <- substr(read_seq,241,320)
  writeLines(c(line1,line2,line3,line4,line5), fileConn)
  close(fileConn)
}

for(read_num in 1:nrow(reads_table)) {
  print(read_num)
  arg_string <- paste("-gapopen 2 -gapextend 1 -penalty -1 -strand minus -query Temp/primer.fasta -subject Temp/r1_",read_num,".fasta -word_size 4 -max_target_seqs 10000000 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \" ", sep="")
  hits_table <- system2("blastn", args = arg_string , stdout = TRUE, stderr = FALSE)
  arg_string2 <- paste("-gapopen 2 -gapextend 1 -penalty -1 -strand minus -query Temp/primer.fasta -subject Temp/r1_",read_num,".fasta -word_size 4 -max_target_seqs 1 -outfmt 0", sep="")
  
  file.remove(paste("rm Temp/r1_",read_num,".fasta",sep="")) # remove files as you go
  
  # to record individual alignments 
  hits_align <- system2("blastn", args = arg_string2, stdout = TRUE, stderr = FALSE)
  # to see alignments as they're generated (can comment out)
  print(hits_align[30])
  print(hits_align[31])
  print(hits_align[32])

  if(length(hits_table) != 0) {
    match <- as.data.frame(t(as.data.frame(strsplit(hits_table, "\t"), stringsAsFactors = FALSE)),row.names = read_num, stringsAsFactors=FALSE)
    colnames(match) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")
    # if more than one alignment want only the first listed alignment (which will be the longest with highest score)
    match <- match[1,] # only look at the best hit
    match[3:14] <- as.numeric(match[3:14]) # convert the number outputs from char to num

    # record all values in table
    R1_insert_trim[read_num,1:2] <- unlist(match[1,1:2])
    R1_insert_trim[read_num,3:14] <- unlist(as.numeric(match[1,3:14]))
    R1_insert_trim[read_num,16] <- primer_seq
    R1_insert_trim[read_num,17] <- hits_align[30]
    R1_insert_trim[read_num,18] <- hits_align[31]
    R1_insert_trim[read_num,19] <- hits_align[32]
    print(match) # uncomment to see the alignment statistics as they're recorded

  } else {print("no matches")
    R1_insert_trim$sseqid[read_num] <- reads_table$Id[read_num]
    R1_insert_trim$slen[read_num] <- nchar(reads_table$R1_seq[read_num])
    R1_insert_trim$qlen[read_num] <- nchar(primer_seq)
  }
}

system("rm -rf Temp/")
system("mkdir Temp")

save(R1_insert_trim, file = paste("Tables/", sample_prefix, "_R1_insert_trim.RData", sep=""))

#######################################################################################################################

#######################################################################################################################
# find adapter sequences at end of read2
#######################################################################################################################

# retreive adapters and indexes for this sample
R2_adapter_trim <- array(NA,c(nrow(reads_table),19))
R2_adapter_trim <- as.data.frame(R2_adapter_trim, stringsAsFactors = FALSE)
colnames(R2_adapter_trim) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen","trim_coord","adapter_seq","align1","align2","align3")
N_residues <- unlist(gregexpr(pattern ="N", sample_seqs$adapter_seq_i7_r1_UMI_end)[1]) # where are the N residues in the adapter sequence
UMIlen <- max(N_residues) - min(N_residues) + 1 # how long is the UMI according to the adapter sequence
Nseq <- substring("NNNNNNNNNNNNNNNNNNNN",1,nchar(reads_table$R3_seq[1])) # how many Ns are being replaced?
adapter_seq <- sample_seqs$adapter_seq_i7_r1_UMI_end # can be found at the 3' end of read 2 in reverse complement orientation

# quick loop to write a fasta file for every adapter
for(read_num in 1:nrow(reads_table)) {
  UMI <- reads_table$R3_seq[read_num]
  UMI <- substring(UMI,1,UMIlen)                  # need this in case the UMI length shorter than the R3 length
  adapter_seq_UMI <- sub(Nseq, UMI, adapter_seq)      # replace N residues with UMI sequence from read 3
  fileConn <- file(paste("Temp/adapter_",read_num,".fasta",sep=""))
  line1 <- ">"
  line2 <- substr(adapter_seq_UMI,1,80)
  line3 <- substr(adapter_seq_UMI,81,160)
  line4 <- substr(adapter_seq_UMI,161,240)
  line5 <- substr(adapter_seq_UMI,241,320)
  writeLines(c(line1,line2,line3,line4,line5), fileConn)
  close(fileConn)
}

# quick loop to write a fasta file for every read2
for(read_num in 1:nrow(reads_table)) {
  fileConn <- file(paste("Temp/r2_",read_num,".fasta",sep=""))
  line1 <- paste(">",reads_table$Id[read_num],sep="")
  read_seq <- reads_table$R2_seq[read_num]
  line2 <- substr(read_seq,1,80)
  line3 <- substr(read_seq,81,160)
  line4 <- substr(read_seq,161,240)
  line5 <- substr(read_seq,241,320)
  writeLines(c(line1,line2,line3,line4,line5), fileConn)
  close(fileConn)
}

for(read_num in 1:nrow(reads_table)) {
  print(read_num)
  arg_string <- paste("-gapopen 2 -gapextend 1 -penalty -1 -strand minus -query Temp/adapter_",read_num,".fasta -subject Temp/r2_",read_num,".fasta -word_size 4 -max_target_seqs 10000000 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \" ", sep="")
  hits_table <- system2("blastn", args = arg_string , stdout = TRUE, stderr = FALSE)
  arg_string2 <- paste("-gapopen 2 -gapextend 1 -penalty -1 -strand minus -query Temp/adapter_",read_num,".fasta -subject Temp/r2_",read_num,".fasta -word_size 4 -max_target_seqs 1 -outfmt 0", sep="")
  
  file.remove(paste("rm Temp/r2_",read_num,".fasta",sep="")) # remove files as you go
  file.remove(paste("rm Temp/adapter_",read_num,".fasta",sep="")) # remove files as you go
  
  # to record individual alignments 
  hits_align <- system2("blastn", args = arg_string2, stdout = TRUE, stderr = FALSE)
  # to see alignments as they're generated (can comment out)
  print(hits_align[30])
  print(hits_align[31])
  print(hits_align[32])

  if(length(hits_table) != 0) {

    match <- as.data.frame(t(as.data.frame(strsplit(hits_table, "\t"), stringsAsFactors = FALSE)),row.names = read_num, stringsAsFactors=FALSE)
    colnames(match) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")
    match <- match[1,] # only look at the best hit
    match[3:14] <- as.numeric(match[3:14]) # convert the number outputs from char to num

    # record all values in table
    R2_adapter_trim[read_num,1:2] <- unlist(match[1,1:2])
    R2_adapter_trim[read_num,3:14] <- unlist(as.numeric(match[1,3:14]))
    R2_adapter_trim[read_num,16] <- adapter_seq
    R2_adapter_trim[read_num,17] <- hits_align[30]
    R2_adapter_trim[read_num,18] <- hits_align[31]
    R2_adapter_trim[read_num,19] <- hits_align[32]
      
    print(match) # uncomment to see the alignment statistics as they're recorded

  } else {print("no matches")
    R2_adapter_trim$sseqid[read_num] <- reads_table$Id[read_num]
    R2_adapter_trim$slen[read_num] <- nchar(reads_table$R2_seq[read_num])
    R2_adapter_trim$qlen[read_num] <- nchar(adapter_seq)
  }
}

system("rm -rf Temp/")
system("mkdir Temp")

save(R2_adapter_trim, file = paste("Tables/", sample_prefix, "_R2_insert_trim.RData", sep=""))

#######################################################################################################################


#######################################################################################################################
# trim reads and store in table
#######################################################################################################################

# R1 matched sequences are only trimmed if the match fall within 16 bases of the adapter/insert sequence start
# and if the bitscore is greater than 14

R1_matches <- which(R1_insert_trim$qlen - R1_insert_trim$qend < 16 & R1_insert_trim$bitscore > 16)

# these criteria for rejected/kept alignments can be reviewed using the following loop to look at individual alignments
R1_mismatches <- setdiff( 1:nrow(R1_insert_trim), R1_matches)

print("R1 trimming that is kept")
for(i in 1:length(R1_matches)) {
  R1 <- R1_matches[i]
  print(R1_insert_trim$align1[R1])
  print(R1_insert_trim$align2[R1])
  print(R1_insert_trim$align3[R1])
}

print("R1 trimming that is discarded")
for(i in 1:length(R1_mismatches)) {
  R1 <- R1_mismatches[i]
  print(R1_insert_trim$align1[R1])
  print(R1_insert_trim$align2[R1])
  print(R1_insert_trim$align3[R1])
}

# read 1 trim coordinates

R1_insert_trim$trim_coord <- R1_insert_trim$slen # the 5' sequence end is the trimming end unless otherwise specified
R1_insert_trim$trim_coord[R1_matches] <- R1_insert_trim$send[R1_matches]

# R2 matched sequences are only trimmed if the match fall within 16 bases of the adapter/insert sequence start
# and if the bitscore is greater than 14

R2_matches <- which(R2_adapter_trim$qlen - R2_adapter_trim$qend < 16 & R2_adapter_trim$bitscore > 16)
# these criteria for rejected/kept alignments can be reviewed using the following loop to look at individual alignments
R2_mismatches <- setdiff( 1:nrow(R2_adapter_trim), R2_matches)

print("R2 trimming that is kept")
for(i in 1:length(R2_matches)) {
  R2 <- R2_matches[i]
  print(R2_adapter_trim$align1[R2])
  print(R2_adapter_trim$align2[R2])
  print(R2_adapter_trim$align3[R2])
}

print("R2 trimming that is discarded")
for(i in 1:length(R2_mismatches)) {
  R2 <- R2_mismatches[i]
  print(R2_adapter_trim$align1[R2])
  print(R2_adapter_trim$align2[R2])
  print(R2_adapter_trim$align3[R2])
}

# read 2 trim coordinates
R2_adapter_trim$trim_coord <- R2_adapter_trim$slen # the 5' sequence end is the trimming end unless otherwise specified
R2_adapter_trim$trim_coord[R2_matches] <- R2_adapter_trim$send[R2_matches]

save(R1_insert_trim, file = paste("Tables/", sample_prefix, "_R1_insert_trim.RData", sep=""))
save(R2_adapter_trim, file = paste("Tables/", sample_prefix, "_R2_insert_trim.RData", sep=""))

# append the coordinates to reads_table
r1trim <- as.data.frame(cbind(R1_insert_trim$sseqid, R1_insert_trim$trim_coord), stringsAsFactors = FALSE)
colnames(r1trim) <- c("Id","R1trim")
reads_table <- merge(reads_table, r1trim, by.x=1, by.y=1, all.x=TRUE)
reads_table$R1trim <- as.integer(reads_table$R1trim)

r2trim <- as.data.frame(cbind(R2_adapter_trim$sseqid, R2_adapter_trim$trim_coord), stringsAsFactors = FALSE)
colnames(r2trim) <- c("Id","R2trim")
reads_table <- merge(reads_table, r2trim, by.x=1, by.y=1, all.x=TRUE)
reads_table$R2trim <- as.integer(reads_table$R2trim)

rm(r1trim, r2trim)

# add four empty columns to reads_table to record the trimmed sequences and quality values
reads_table <- cbind(reads_table, as.data.frame(array(NA, dim = c(nrow(reads_table), 4)), stringsAsFactors = FALSE))
colnames(reads_table)[15:18] <- c("R1trimseq", "R1trimqual", "R2trimseq", "R2trimqual")

# store trimmed sequences and quality values in reads_table
for (read_num in 1:nrow(reads_table)) {

  # if any read has no hits for trimming, use untrimmed end values
  if (is.na(reads_table$R1trim[read_num])) {
    reads_table$R1trim[read_num] <- nchar(as.character(reads_table$R1_seq[read_num]))
  }
  if (is.na(reads_table$R2trim[read_num])) {
    reads_table$R2trim[read_num] <- nchar(as.character(reads_table$R2_seq[read_num]))
  }

  # store trim values in reads_table
  reads_table$R1trimseq[read_num] <- substr(reads_table$R1_seq[read_num], start = 1 , stop = reads_table$R1trim[read_num])
  reads_table$R1trimqual[read_num] <- substr(reads_table$R1_qual[read_num], start = 1 , stop = reads_table$R1trim[read_num])
  reads_table$R2trimseq[read_num] <- substr(reads_table$R2_seq[read_num], start = reads_table$R2trim_5p[read_num] , stop = reads_table$R2trim[read_num])
  reads_table$R2trimqual[read_num] <- substr(reads_table$R2_qual[read_num], start = reads_table$R2trim_5p[read_num] , stop = reads_table$R2trim[read_num])
}

# make trimmed version of read1 & read2 shortread variables

read1trim <- read1[1:nrow(reads_table)] # make a copy of read1 with correct number of rows (content doesn't matter)
read2trim <- read2[1:nrow(reads_table)] # make a copy of read2 with correct number of rows (content doesn't matter)
read1trim[1:nrow(reads_table)] <- ShortReadQ( sread = DNAStringSet(""), quality = BStringSet(""), id = BStringSet("")) # empty all reads
read2trim[1:nrow(reads_table)] <- ShortReadQ( sread = DNAStringSet(""), quality = BStringSet(""), id = BStringSet("")) # empty all reads

r1_id_long = BStringSet(reads_table$R1_ID_long)
r1_trimseq <- DNAStringSet(reads_table$R1trimseq)
r1_trimquality  <- BStringSet(reads_table$R1trimqual)

r2_id_long = BStringSet(reads_table$R2_ID_long)
r2_trimseq <- DNAStringSet(reads_table$R2trimseq)
r2_trimquality  <- BStringSet(reads_table$R2trimqual)

read1trim <- ShortReadQ(sread = r1_trimseq, quality = r1_trimquality, id = r1_id_long)
read2trim <- ShortReadQ(sread = r2_trimseq, quality = r2_trimquality, id = r2_id_long)

# write fastq files for trimmed reads 
# writeFastq won't overwrite existing files so place files in Temp directory then move to overwrite
writeFastq(read1trim, paste("Temp/", sample_prefix, "_R1trim.fastq", sep="") , mode="w", compress = FALSE)
writeFastq(read2trim, paste("Temp/", sample_prefix, "_R2trim.fastq", sep="") , mode="w", compress = FALSE)
system("mv Temp/*trim.fastq Trimmed")

save(reads_table, file = paste("Tables/",sample_prefix,"_reads_table_trim.RData",sep=""))

# optional: uncomment to see untrimmed and trimmed reads lined up with each other
for(read_num in 1:nrow(reads_table)) {
  print(read_num)
  id_long <- r1_id_long[read_num]
  print(paste(id_long,"read1"))
  print(as.character(sread( read1[which(id(read1) == id_long)] )))
  print(as.character(sread( read1trim[which(id(read1trim) == id_long)] )))
  print(paste(id_long,"read2"))
  id_long <- r2_id_long[read_num]
  print(as.character(sread( read2[which(id(read2) == id_long)] )))
  print(paste(strrep(" ",reads_table$R2trim_5p[read_num]-1), as.character(sread( read2trim[which(id(read2trim) == id_long)] )),sep="") )
}
