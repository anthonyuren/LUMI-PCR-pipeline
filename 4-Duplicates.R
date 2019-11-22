#######################################################################################################################
# LUMI-PCR Script for custom trimming of reads
#######################################################################################################################
args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1] 
# e.g. work_dir <- '/Volumes/ILLUMINA/M-1-24b/'
# e.g. work_dir <- '~/Desktop/M-1-24c/'
setwd(work_dir)

# load a table listing details of the samples that each library came from
sample_table <- read.delim("sample_table.txt", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)

#######################################################################################################################
# load all the insert files and merge them into a single table
# first option is table without UMI sequences
# second option is the full table which can strain memory for tables with millions of inserts
# if UMI sequences will not be analysed for contamination analyses this can be ignored and the short table can be used
#######################################################################################################################
source("4-load_all_inserts_short.R") # creates short table all_inserts
source("4-load_all_inserts_long.R")  # creates long table all_inserts_full, comment out as necessary

# replace all_inserts with all_inserts_full for subsequent steps, comment out as necessary
all_inserts <- all_inserts_full
rm(all_inserts_full)

#######################################################################################################################
# identify duplicate inserts from different samples
#######################################################################################################################

# this inserts are batched by chromosomes and by strand and then clustered to identify duplicates
cluster_width <- 1 # a larger value is less stringent in calling adjacent inserts as duplicates

for(chr_strand in unique(all_inserts$chr_strand)) {
  insert_rows <- which(all_inserts$chr_strand == chr_strand)
  if(length(insert_rows) > 1) {
    dist_mat <-dist(all_inserts$average_base[insert_rows], method="euclidian")
    hclust_inserts <- hclust(dist_mat, method="single")    # hierarchical clustering of the distance matrix
    duplicates <-cutree(hclust_inserts,h = cluster_width)
  } else { duplicates <- 1 }
  all_inserts$duplicates[insert_rows] <- paste(all_inserts$chr_strand[insert_rows], duplicates, sep="_")
  all_inserts$duplicate_coords[insert_rows] <- paste(all_inserts$chr[insert_rows], ":", all_inserts$average_base[insert_rows], all_inserts$strand[insert_rows], "_", duplicates, sep="")
}

# delete the duplicate ids for inserts that are unique
for(i in 1:nrow(all_inserts)) {
  if( length(which(all_inserts$duplicates == all_inserts$duplicates[i] ) ) ==1 ) { 
    all_inserts$duplicates[i] <- ""
    all_inserts$duplicate_coords[i] <- ""}
}

write.table(all_inserts, file = "Tables/all_inserts.txt", quote=FALSE, row.names=FALSE, sep="\t")

#######################################################################################################################
# make a matrix of the duplicate groups vs samples for analysis/visualization
#######################################################################################################################

dup_inserts    <- all_inserts[(all_inserts$duplicates != ""),] # these may require filtering
unique_inserts <- all_inserts[(all_inserts$duplicates == ""),] # these are finished product and don't need further filtering

samples <- unique(dup_inserts$sample)
dup_groups <- unique(dup_inserts$duplicates)

# create a matrix to store the CLONALITY values of each of the duplicate values in each sample
duplicate_array <- as.data.frame(array(data = 0, dim = c(length(samples), length(dup_groups))))
rownames(duplicate_array) <- sort(samples)
colnames(duplicate_array) <- sort(dup_groups)

for(dup in 1:length(dup_groups)) {
  this_dup <- dup_groups[dup]
  for(sam in 1:length(samples)) {
    this_sam <- samples[sam]
    cell <- as.numeric(subset(dup_inserts, dup_inserts$sample == this_sam & dup_inserts$duplicates == this_dup)[12]) # 12 is CLONALITY not normalized clonality
    print(paste(this_dup, this_sam,cell))
    if(is.na(cell)) { cell <-0}
    duplicate_array[this_sam, this_dup] <- cell
  }
}

# list the coordinates for each of the dup groups so they can be included in the dup_group name
dup_groups_bases <- c()
dup_groups_names <- data.frame(dup_groups = rep("",length(dup_groups)) ,
                               bases = rep("",length(dup_groups)) ,
                               names = rep("",length(dup_groups)) ,
                               stringsAsFactors = FALSE)

for(i in 1:length(dup_groups)) {
  dup_groups_names$dup_groups[i] <- dup_groups[i]
  dup_groups_names$bases[i] <- paste( c(unique(dup_inserts$average_base[which(dup_inserts$duplicates == dup_groups[i])])) , collapse = "," )
  name <- dup_groups[i]
  name <- sub(":", paste(":",dup_groups_names$bases[i],"_",sep=""), name)
  dup_groups_names$names[i] <- name
  # rename the colums with the longer string
  replace <- which(colnames(duplicate_array)==dup_groups[i])
  colnames(duplicate_array)[replace] <- name
}

# can be useful to add sort rows/columns by average/min/max values

#######################################################################################################################
# flag duplicates according to whether they are expected
#######################################################################################################################

sample_table <- read.table("sample_table.txt", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)

sample_table$replicate_group
sample_table$overlap_groups

# otherwise they may represent cross contamination/mispriming/PCR artefacts

# list samples that are controls (i.e. ones with only a zero in the overlap column)
controls <- sample_table$sample[which(sample_table$overlap_groups == "0")]

# which regions of the genome are flagged for rescue from duplicate filtering because they are known to be highly recurrent between samples
# e.g. in MuLV infected lymphomas the 3'UTR of Mycn has many sense orientation inserts (sometimes multiple single read inserts per tumor sample)
# additional lines can be added by rbind as needed
rescue_loci <- c("chr12", 12936980, 12937000)
rescue_loci <- rbind(rescue_loci, c("chr12", 12936960, 12936970))

for(dup in dup_groups) { # for each duplicate group
  print(dup)
  dup_samples <- which(dup_inserts$duplicates == dup)
  print(dup_inserts$sample[dup_samples])
  
  # is this duplicate_group ever observed in control samples i.e. overlap_groups == "0"
  if(length(intersect(controls,dup_inserts$sample[dup_samples])) > 1) { 
    print("found in control samples")
    print(intersect(controls, dup_inserts$sample[dup_samples]))
    dup_inserts$control[dup_samples] <- paste(length(intersect(controls, dup_inserts$sample[dup_samples]))) # record number of control samples for all inserts
  }
  
  # which samples are replicates and if these duplicates are limited to a single replicate_group they can be ignored
  rep_list <- sample_table$replicate_group[which(sample_table$sample %in% dup_inserts$sample[dup_samples])]
  print(rep_list)
  dup_inserts$replicate_groups[dup_samples] <- toString(unique(rep_list) ) # record replicate_groups this duplicate group can be found in
  
  # which duplicates are expected becaue they result from mixed DNAs or arise from the same tissue or same organism
  # i.e. which duplicates share overlap_groups
  overlap_list <- sample_table$overlap_groups[which(sample_table$sample %in% dup_inserts$sample[dup_samples])]
  print(overlap_list)
  dup_inserts$expected[dup_samples] <- toString(unique(overlap_list) ) #record overlap_groups this duplicate group can be found in

  # which duplicates fall within regions that should be flagged for rescue from filtering
  for (rrow in 1:nrow(rescue_loci)) {
    rescue_chr <- rescue_loci[rrow,1]
    rescue_start <- rescue_loci[rrow,2]
    rescue_end <- rescue_loci[rrow,3]
    dup_chr <- unique(dup_inserts$chr[dup_samples])
    dup_base <- unique(dup_inserts$average_base[dup_samples])
    if (rescue_chr == dup_chr){
      if (dup_base >= rescue_start & dup_base <= rescue_end) {
             print(paste(rescue_chr,":",rescue_start,"-",rescue_end,sep=""))
        dup_inserts$rescue[dup_samples] <- paste(rescue_chr,":",rescue_start,"-",rescue_end,sep="")
      }
    }
  }
}

flagged_inserts <- rbind(dup_inserts, unique_inserts)
save(flagged_inserts, file = "Tables/flagged_inserts.RData")
write.table(flagged_inserts, file = "Tables/flagged_inserts.txt", quote=FALSE, row.names=FALSE, sep="\t")

  
#########################################################################################################################
# calculate entropy from the first 50 most clonal inserts
#########################################################################################################################

for (this_sample in samples) {
  inserts <- subset(all_inserts, sample == this_sample)
  inserts <- inserts[order(inserts$normclon_collUMIs, decreasing=TRUE),]
  clonalities <- inserts$normclon_collUMIs
  print(this_sample)
  
  ent <- 0.0
  unique(clonalities[1:50])
  for (clon in clonalities[1:50]) {
    pi <- clon/sum(clonalities[1:50])
    ent <- ent + (pi * log(1/pi))
  }
  print(ent)
  clon_freqs <- clonalities[1:50]/sum(clonalities[1:50])
  print(-sum(clon_freqs * log(clon_freqs)))
  
  ent <- round(ent,digits = 3)
  barplot(clonalities[1:50], main= paste(as.character(this_sample)," E=",as.character(ent),sep=""),
          xlab="Normalized Clonality", col="Red")
}

