#######################################################################################################################
# LUMI-PCR Script for samtools and mapping
#######################################################################################################################

library(ShortRead)
library(Biostrings)
library(GenomicAlignments)
sessionInfo()

#######################################################################################################################
# Input filenames/arguments
#######################################################################################################################

args <- commandArgs(trailingOnly = TRUE)
sample_prefix <- args[1] # e.g. sample_prefix <- "M-1_S95" 
                         # e.g. sample_prefix <- "L-6_S78"
work_dir <- args[2]      # e.g. work_dir <- '/Volumes/ILLUMINA/M-1-24/'
                         # e.g. work_dir <- '/Volumes/ILLUMINA/L-5-10/'
genome <- args[3]        # e.g. genome <- "mm10"
                         # e.g. genome <- "hg38"
insert_seq <- args[4]    # e.g. insert_seq <- "MMuLV"
                         # e.g. insert_seq <- "Lenti"

setwd(work_dir)

R1file <- paste("Trimmed/", sample_prefix, "_R1trim.fastq", sep="")
R2file <- paste("Trimmed/", sample_prefix, "_R2trim.fastq", sep="")
table_file <- paste("Tables/", sample_prefix, "_reads_table_trim.RData", sep="")

print(work_dir)
print("Read files.")
print(R1file)
print(R2file)
print(table_file)


#######################################################################################################################
# map reads to genome and run samtools
#######################################################################################################################
arg_string <- paste("-query ", R1file, " -query_mate ", R2file, " -db ../Genomes/", genome, " -paired -infmt fastq -outfmt sam -splice F > Mapped/", sample_prefix, ".sam", sep = "")
system2("magicblast", args = arg_string, stdout = TRUE, stderr = TRUE)

system(paste("samtools view -f 0x2 Mapped/", sample_prefix, ".sam > Mapped/", sample_prefix, "_paired.sam", sep="")) # keep only the properly paired mappings and store in sam file
sam_file <- readLines( paste("Mapped/", sample_prefix, ".sam", sep="") )
pair_sam_file <- readLines( paste("Mapped/", sample_prefix, "_paired.sam", sep="") )
head_sam_file <- sam_file[grep("@", sam_file)]                                        # get the sam_file header to add to the paired_sam_file

chrsam_file <- pair_sam_file[grep("=", pair_sam_file)]                                # only keep pairs that map to the same chromosomes (others can't be interpreted)
writeLines(c(head_sam_file, chrsam_file), paste("Mapped/", sample_prefix, "_paired.sam", sep="") )
system(paste("samtools view -Sb Mapped/", sample_prefix, "_paired.sam > Mapped/", sample_prefix, "_paired.bam", sep=""))    # convert to bam file
system(paste("samtools sort Mapped/", sample_prefix, "_paired.bam -o Mapped/", sample_prefix, "_paired_sorted.bam", sep=""))# sort bam file
system(paste("samtools index Mapped/", sample_prefix, "_paired_sorted.bam", sep=""))                                        # create index for sorted bam file


#######################################################################################################################
# map reads to insert sequence and run samtools
#######################################################################################################################
arg_string <- paste("-query ", R1file, " -query_mate ", R2file, " -db ../Genomes/", insert_seq, " -paired -infmt fastq -outfmt sam -splice F > Mapped/", sample_prefix, "_",insert_seq, ".sam", sep = "")

system2("magicblast", args = arg_string, stdout = TRUE, stderr = TRUE)

# comparing to mapping of untrimmed reads to trimmed reads gives good idea of the effectiveness of trimming
# more reads should map post trimming and there should be less soft clipping in the cigar string

insert_prefix <- paste(sample_prefix,"_",insert_seq, sep="")

system(paste("samtools view -f 0x2 Mapped/", insert_prefix, ".sam > Mapped/", insert_prefix, "_paired.sam", sep=""))      # keep only the properly paired mappings and store in sam file
sam_file <- readLines( paste("Mapped/", insert_prefix, ".sam", sep="") )
pair_sam_file <- readLines( paste("Mapped/", insert_prefix, "_paired.sam", sep="") )
head_sam_file <- head_sam_file  # keep this the same
chrsam_file <- pair_sam_file

writeLines(c(head_sam_file, chrsam_file), paste("Mapped/", insert_prefix, "_paired.sam", sep="") )

system(paste("samtools view -Sb Mapped/", insert_prefix, "_paired.sam > Mapped/", insert_prefix, "_paired.bam", sep=""))    # convert to bam file
system(paste("samtools sort Mapped/", insert_prefix, "_paired.bam -o Mapped/", insert_prefix, "_paired_sorted.bam", sep=""))# sort bam file
system(paste("samtools index Mapped/", insert_prefix, "_paired_sorted.bam", sep=""))                                        # create index for sorted bam file


#######################################################################################################################
# import bamfiles as table
#######################################################################################################################

# create GAlignmentPairs object from the sorted bam file
genome_align <- readGAlignmentPairs(paste("Mapped/", sample_prefix, "_paired_sorted.bam", sep=""), index=paste("Mapped/", sample_prefix, "_paired_sorted.bam",sep=""), use.names=TRUE, param=ScanBamParam(what=scanBamWhat()), strandMode = 1)


# qstart is always first relative to qend for both reads
# i.e. start end coordinates read right to left on genome and orientation is determined by strand
# an integration of MuLV on the forward strand using primers that amplify the 5' insert-genome junction
# should have a "+" strand for read 1 and a "-" strand for read 2

genome_map <- data.frame(Id=character(),
                         chr=character(),
                         strandr1=character(),
                         strandr2=character(),
                         widthr1=integer(),
                         widthr2=integer(),
                         qwidthr1=integer(),
                         qwidthr2=integer(),
                         qstartr1=integer(),
                         qstartr2=integer(),
                         qendr1=integer(),
                         qendr2=integer(),
                         CIGARr1=character(),
                         CIGARr2=character(),
                         isizer1=integer(),
                         isizer2=integer(),
                         start=integer(),
                         end=integer(), stringsAsFactors = FALSE)

if(length(genome_align[]) >0 ) {
for (pair_num in 1:length(genome_align)) {
  genome_map[pair_num,] <- rep(NA,18)
  pair_df <- elementMetadata(genome_align[[pair_num]])
  genome_map$Id[pair_num]       <- pair_df$qname[1]
  genome_map$chr[pair_num]      <- as.character(pair_df$rname[1])
  genome_map$strandr1[pair_num] <- as.character(pair_df$strand[1])
  genome_map$strandr2[pair_num] <- as.character(pair_df$strand[2])
  genome_map$widthr1[pair_num]  <- width(genome_align[[pair_num]])[1]
  genome_map$widthr2[pair_num]  <- width(genome_align[[pair_num]])[2]
  genome_map$qwidthr1[pair_num] <- pair_df$qwidth[1]
  genome_map$qwidthr2[pair_num] <- pair_df$qwidth[2]
  genome_map$qstartr1[pair_num] <- start(genome_align[[pair_num]])[1]
  genome_map$qstartr2[pair_num] <- start(genome_align[[pair_num]])[2]
  genome_map$qendr1[pair_num]   <- end(genome_align[[pair_num]])[1]
  genome_map$qendr2[pair_num]   <- end(genome_align[[pair_num]])[2]
  genome_map$CIGARr1[pair_num]  <- pair_df$cigar[1]
  genome_map$CIGARr2[pair_num]  <- pair_df$cigar[2]
  genome_map$isizer1[pair_num]  <- pair_df$isize[1]
  genome_map$isizer2[pair_num]  <- pair_df$isize[2]
  genome_map$start[pair_num]    <- min(genome_map$qstartr1[pair_num], genome_map$qstartr2[pair_num])
  genome_map$end[pair_num]      <- max(genome_map$qendr1[pair_num], genome_map$qendr2[pair_num])
}
}

insert_align <- readGAlignmentPairs(paste("Mapped/", sample_prefix,"_",insert_seq, "_paired_sorted.bam",sep=""), index=paste("Mapped/", sample_prefix,"_",insert_seq, "_paired_sorted.bam",sep=""), use.names=TRUE, param=ScanBamParam(what=scanBamWhat()), strandMode = 1)

# convert GAlignmentPairs object into a table, the i stands for insert mapping
insert_map <- data.frame(Id=character(),
                         chri=character(),
                         strandr1i=character(),
                         strandr2i=character(),
                         widthr1i=integer(),
                         widthr2i=integer(),
                         qwidthr1i=integer(),
                         qwidthr2i=integer(),
                         qstartr1i=integer(),
                         qstartr2i=integer(),
                         qendr1i=integer(),
                         qendr2i=integer(),
                         CIGARr1i=character(),
                         CIGARr2i=character(),
                         isizer1i=integer(),
                         isizer2i=integer(),
                         starti=integer(),
                         endi=integer(), stringsAsFactors = FALSE)

if(length(insert_align[]) >0 ) {
for (pair_num in 1:length(insert_align[])) {
  insert_map[pair_num,] <- rep(NA,18)
  pair_df <- elementMetadata(insert_align[[pair_num]])
  insert_map$Id[pair_num]        <- pair_df$qname[1]
  insert_map$chri[pair_num]      <- as.character(pair_df$rname[1])
  insert_map$strandr1i[pair_num] <- as.character(pair_df$strand[1])
  insert_map$strandr2i[pair_num] <- as.character(pair_df$strand[2])
  insert_map$widthr1i[pair_num]  <- width(insert_align[[pair_num]])[1]
  insert_map$widthr2i[pair_num]  <- width(insert_align[[pair_num]])[2]
  insert_map$qwidthr1i[pair_num] <- pair_df$qwidth[1]
  insert_map$qwidthr2i[pair_num] <- pair_df$qwidth[2]
  insert_map$qstartr1i[pair_num] <- start(insert_align[[pair_num]])[1]
  insert_map$qstartr2i[pair_num] <- start(insert_align[[pair_num]])[2]
  insert_map$qendr1i[pair_num]   <- end(insert_align[[pair_num]])[1]
  insert_map$qendr2i[pair_num]   <- end(insert_align[[pair_num]])[2]
  insert_map$CIGARr1i[pair_num]  <- pair_df$cigar[1]
  insert_map$CIGARr2i[pair_num]  <- pair_df$cigar[2]
  insert_map$isizer1i[pair_num]  <- pair_df$isize[1]
  insert_map$isizer2i[pair_num]  <- pair_df$isize[2]
  insert_map$starti[pair_num]    <- min(insert_map$qstartr1[pair_num], insert_map$qstartr2[pair_num])
  insert_map$endi[pair_num]      <- max(insert_map$qendr1[pair_num], insert_map$qendr2[pair_num])
}
}
# save both mapping tables
save(insert_map, genome_map, file=paste("Tables/", sample_prefix, "_mapping.RData", sep=""))

#######################################################################################################################
# add columns to reads_table
#######################################################################################################################

# remove rows that have improbably sized mapping i.e. > 1000bp separation between reads
insert_map <- insert_map[which(abs(insert_map$isizer1i) < 1000),]
genome_map <- genome_map[which(abs(genome_map$isizer1) < 1000),]

load(table_file)
reads_table <- merge(reads_table, insert_map[,c(1,2,3,15,16,17,18)], by.x=1, by.y=1, all.x=TRUE)
reads_table <- merge(reads_table, genome_map[,c(1,2,3,15,16,17,18)], by.x=1, by.y=1, all.x=TRUE)

reads_table <- cbind(reads_table, NA)
colnames(reads_table)[31] <- "genome_or_insert"
reads_table$genome_or_insert[which(!is.na(reads_table$chri))] <- "insert"
reads_table$genome_or_insert[which(!is.na(reads_table$chr))]  <- "genome"

# if the read pair maps to both insert and genome (this is rare if it happens at all)
# call both and label which has the longer mapping or default to genome if lengths are same

map_both <- intersect(which(!is.na(reads_table$isizer1)), which(!is.na(reads_table$isizer1i)))

for (read in map_both) {
  if (abs(reads_table$isizer1[read]) < abs(reads_table$isizer1i[read])) {reads_table$genome_or_insert[read] <- "both_insert"}
  else {reads_table$genome_or_insert[read] <- "both_genome"}
}

reads_table <- cbind(reads_table, NA)
colnames(reads_table)[28] <- "insert_end"
reads_table$insert_end[which(reads_table$strandr1 == "+")] <- reads_table$end[which(reads_table$strandr1 == "+")]
reads_table$insert_end[which(reads_table$strandr1 == "-")] <- reads_table$start[which(reads_table$strandr1 == "-")]

save(reads_table, file=paste("Tables/", sample_prefix, "_mapped_reads_table.RData", sep=""))

#######################################################################################################################
