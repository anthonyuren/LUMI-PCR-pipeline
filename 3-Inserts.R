#######################################################################################################################
# LUMI-PCR Script for grouping reads that map together into inserts and counting fragments/UMIs/reads per insert
#######################################################################################################################

library("ShortRead")
sessionInfo()

args <- commandArgs(trailingOnly = TRUE)
sample_prefix <- args[1] 
work_dir <- args[2]
# e.g. sample_prefix <- "M-1_S95" ;  work_dir <- '/Volumes/ILLUMINA/M-1-24/'
# e.g. sample_prefix <- "M-13_S107" ;  work_dir <- '/Volumes/ILLUMINA/M-1-24/'
# e.g. sample_prefix <- "M-21_S115" ;  work_dir <- '/Volumes/ILLUMINA/M-1-24/'

setwd(work_dir)
print(work_dir)

#######################################################################################################################
# Select only rows that map to genome and have expected insert-genome junction
#######################################################################################################################

table_file <- paste("Tables/", sample_prefix, "_mapped_reads_table.RData", sep="")
load(table_file)
good_reads_table <- reads_table
good_reads_table <- good_reads_table[union(which(good_reads_table$genome_or_insert == "genome"), which(good_reads_table$genome_or_insert == "both_genome")), ]
# could optionally include other selective criteria here e.g. read quality, UMI quality, read length, mapping length, 
save(good_reads_table, file = paste("Tables/", sample_prefix, "_good_reads_table.RData", sep=""))
# to group reads into inserts will only need a subset of fields 

grouped_reads <- subset(good_reads_table, select = c("Id","chr","insert_end","strandr1","start","end","R3_seq","R2_seq"))
rm(table_file, reads_table, good_reads_table)

#######################################################################################################################
# Assign each insert to a group based on chromosome and strand
#######################################################################################################################

chr_strand_groups <- paste(grouped_reads$chr, grouped_reads$strandr1, sep="")
grouped_reads <- cbind(grouped_reads, chr_strand_groups)
grouped_reads <- grouped_reads[order(grouped_reads$chr, grouped_reads$strand, grouped_reads$insert_end),] # make pretty
rm(chr_strand_groups)

# add a column for insert groups
grouped_reads <- cbind(grouped_reads, NA)
colnames(grouped_reads)[10] <- "insert_group"

#######################################################################################################################
# For each chromosome/strand group of inserts, cluster insertions within 10bp of each other
#######################################################################################################################

chr_strand_list <- unique(grouped_reads$chr_strand_groups) # lists all the chromosome-strand combinations observed

for (chr_str in chr_strand_list) { # e.g chr_str <- "chr19+"
  # reducing list of insert positions to unique positions cuts down processing/memory required for dist matrix & clustering
  insert_positions <- unique(grouped_reads$insert_end[grouped_reads$chr_strand_groups == chr_str])
  print(paste(chr_str,insert_positions,sep=":"))
  dist_mat <- NA
  dist_mat <- dist(insert_positions, method = "euclidean") # if there's only one insert_position this is empty i.e. length == 0
  
  if(length(dist_mat) > 0 ) {                              # if there's more than one insert position
    hclust_inserts <- hclust(dist_mat, method="single")    # hierarchical clustering of the distance matrix
    groups <-cutree(hclust_inserts,h = 10)                 # cut the cluster tree with 10bp window
  } else {                                                 # else there is only one insert position
    groups <- 1                                            # all reads at that position belong to group 1
  }
  
  groups <- paste(chr_str, groups, sep="_")
  group_list <- cbind(insert_positions, groups)
  colnames(group_list) <- c("insert_positions", "groups")
  group_list <- as.data.frame(group_list, stringsAsFactors = FALSE)
  print(group_list) # uncomment to monitor grouping
  
  for(insert_pos in group_list$insert_positions) {          # e.g. insert_pos <- 45476625
    this_group <- group_list$groups[which(group_list$insert_positions == insert_pos)]
    grouped_reads$insert_group[which(grouped_reads$chr_strand_groups == chr_str & grouped_reads$insert_end == insert_pos)] <- this_group
  }
}


# screen out recombinant UMIs i.e. identical UMIs between different inserts
grouped_reads <- cbind(grouped_reads, rep(0,nrow(grouped_reads)))
colnames(grouped_reads)[11] <- "UMI_flag"

for (UMI in unique(grouped_reads$R3_seq)) {
  UMIreads <- grouped_reads[which(grouped_reads$R3_seq == UMI),]
  # check for recombinant UMIs i.e. UMIs that are the same between different integrations
  # print(paste(UMI, nrow(UMIreads), unique(UMIreads$insert_group)))
  # these all contain many G residues
  if(length(unique(UMIreads$insert_group)) > 2) {grouped_reads[which(grouped_reads$R3_seq == UMI), 11] <- 1} #flag these for removal
  
  # check for identical UMIs with differing fragment lengths
  # print(paste(UMI, nrow(UMIreads), unique(UMIreads$end - UMIreads$start)))
  # nearly all have the same UMI between different inserts, 
  # very occasionally fragments from same insert with identical UMI will differ by up to 5bp and these are not removed
}

print(paste(sum(grouped_reads$UMI_flag) , "/", nrow(grouped_reads), "reads eliminated due to UMI crossover =", 100*sum(grouped_reads$UMI_flag)/nrow(grouped_reads), "%" ))
grouped_reads <- grouped_reads[which(grouped_reads$UMI_flag == 0) , 1:10]
rm(hclust_inserts, dist_mat, groups, group_list, this_group, insert_positions, insert_pos, chr_str, chr_strand_list)

save(grouped_reads, file= paste("Tables/", sample_prefix, "_grouped_reads.RData", sep=""))

#######################################################################################################################
# function to collapse UMIs that differ by 1bp 
#######################################################################################################################

collapse_function <- function(allUMIs) {
  # test allUMIs <- grouped_reads$R3_seq 
  uniqueUMIs <- as.data.frame(unique(allUMIs), stringsAsFactors = FALSE)
  uniqueUMIs <- cbind(uniqueUMIs, rep(0, length(uniqueUMIs)))
  colnames(uniqueUMIs) <- c("UMI", "count")
  
  for (elem in 1:nrow(uniqueUMIs)) {
    UMI <- uniqueUMIs[elem,1]
    uniqueUMIs[elem,2] <- length(which(allUMIs == UMI)) #how many times does this unique UMI appear in the full allUMIs list
  }
  
  # sort list so the most abundant one is first
  currentUMIs <- as.data.frame( cbind(uniqueUMIs[order(uniqueUMIs$count,decreasing = TRUE) ,], rep(0,nrow(uniqueUMIs))) , stringsAsFactors = FALSE) 
  colnames(currentUMIs) <- c("UMI", "count", "collapsed_count")
  
  for (elem in 1:nrow(currentUMIs) ) { # for each UMI in the list
    UMI <- currentUMIs[elem,1]
    # print(UMI)
    coords <- which(adist(UMI, currentUMIs[,1]) < 2) # equivalent to Hamming distance 1
    
    if (currentUMIs[elem,2] > 0){ # if the counts of this UMI row have not already been reassigned to another UMI
      if (length(coords) > 1) {
        currentUMIs[coords[1],3] <- sum(as.numeric( currentUMIs[coords[1],2] ), as.numeric( currentUMIs[coords[-1],2] ) ) #sum merged UMIs in column 3
        currentUMIs[coords[1],2] <- 0 #empty column 2
        currentUMIs[coords[-1],] <- 0 #set all columns merged UMI sequences to zero
      } else {
        currentUMIs[coords[1],3] <- currentUMIs[coords[1],2]
        currentUMIs[coords[1],2] <- 0
      }
    }
    # print(currentUMIs)
  }
  return(currentUMIs)
}

# test function on all the UMIs in the sample
# collapse_all_UMIs_in_sample <- collapse_function(grouped_reads$R3_seq) # test function on full list of UMIs
# length(unique(grouped_reads$R3_seq))
# length(which(collapse_all_UMIs_in_sample[,1] != 0))

# test function on an individual set of grouped reads
# collapse_function(grouped_reads$R3_seq[which(grouped_reads$insert_group == "chr12-_1")])

#######################################################################################################################
# Make the insert table
#######################################################################################################################

insert_table <- data.frame(sample=character(),               # sample id
                           insert_group=character(),         # unique group id for reads mapping to this insert
                           insert_string=character(),        # unique insert string
                           chr=character(),                  # chromosome
                           insert_base=numeric(),            # insert-genome junction (best estimate if mapping ambiguous reads)
                           strand=character(),               # strand/orientation
                           fragments=integer(),              # count of sheared fragment ends
                           collapsed_UMIs=integer(),         # count of UMIs after collapsing with hamming distance = 1
                           UMIs=integer(),                   # count of unique UMIs
                           reads=integer(),                  # count of reads
                           clon_fragments=numeric(),         # clonality values of the counts
                           clon_collUMIs=numeric(),
                           clon_UMIs=numeric(),
                           clon_reads=numeric(),
                           normclon_fragments=numeric(),     # normalized clonality values of the counts
                           normclon_collUMIs=numeric(),
                           normclon_UMIs=numeric(),
                           normclon_reads=numeric(),
                           insert_positions=character(),     # all insert-genome junction positions for these reads
                           first_pos=integer(),              # junction closest to insert end
                           last_pos=integer(),               # junction furthest from insert end
                           frag_counts=character(),          # count of fragments for each position
                           coll_UMI_counts=character(),      # count of collapsed UMIs for each position
                           read_counts=character(),          # count of reads for each position
                           UMIseqs=character(), stringsAsFactors = FALSE)

insert_list <- unique(grouped_reads$insert_group)
length(insert_list)

for (i in 1:length(insert_list)) {
  
  insert_table[i,] <- rep(NA,25)
  insert <- insert_list[i]
  insert_reads <- grouped_reads[which(grouped_reads$insert_group == insert),]

  insert_table$sample[i]            <- sample_prefix
  insert_table$insert_group[i]      <- insert_list[i]
  insert_table$chr[i]               <- unique(insert_reads$chr)
  insert_table$strand[i]            <- unique(insert_reads$strandr1)
  
  # for inserts that grouop reads mapping to multiple insert-genome junctions
  # record the counts for each junction
  insert_table$insert_positions[i]  <- toString(unique(insert_reads$insert_end))
  if (insert_table$strand[i] == "+") {                                    # on the for/+ strand
  insert_table$fragments[i]         <- length(unique(insert_reads$start)) # $start is the sheared end
  insert_table$first_pos[i]         <- min(insert_reads$end)              # $end is the insert end
  insert_table$last_pos[i]          <- max(insert_reads$end)
  } else {                                                                # on the rev/- strand 
    insert_table$fragments[i]       <- length(unique(insert_reads$end))   # $end is the sheared end
    insert_table$first_pos[i]       <- min(insert_reads$start)            # $start is the insert end
    insert_table$last_pos[i]        <- max(insert_reads$start)
  }
  insert_table$UMIs[i]              <- length(unique(insert_reads$R3_seq))
  collapsed_UMIs                    <- collapse_function(insert_reads$R3_seq)
  insert_table$collapsed_UMIs[i]    <- length(which(collapsed_UMIs[,3] >0))
  insert_table$reads[i]             <- nrow(insert_reads)

  # record remaining fields as strings containing many rows of input
  # if an insert has thousands of reads these strings can get very, very long
  # can comment out UMIseqs field if necessary
  UMIs <- list()
  UMI_count <- list()
  read_count <- list()
  frag_count <- list()
  for (ipos in 1:length( unique(insert_reads$insert_end) )) {
    this_pos <- unique(insert_reads$insert_end)[ipos]
    this_pos_reads <- insert_reads[which(insert_reads$insert_end == this_pos),]
    
    if (insert_table$strand[i] == "+") { # count sheared fragments for each position
      frag_count[ipos] <- as.character(length(unique(this_pos_reads$start)))
    } else {
      frag_count[ipos] <- as.character(length(unique(this_pos_reads$end)))
    }
    read_count[ipos] <- as.character(nrow(this_pos_reads))
    collapsed_UMIs   <- collapse_function(this_pos_reads$R3_seq)
    UMI_count[ipos]  <- as.character(length(which(collapsed_UMIs[,3] >0)))
    UMIs[ipos]       <- toString(this_pos_reads$R3_seq)
  }

  insert_table$frag_counts[i]      <- paste(frag_count, sep="", collapse = ",")
  insert_table$coll_UMI_counts[i]  <- paste(UMI_count, sep="", collapse = ",")
  insert_table$read_counts[i]      <- paste(read_count, sep="", collapse =",")
  # print(UMIs)
  # might want to comment out UMI list if there are too many reads
  insert_table$UMIseqs[i]          <- paste(UMIs ,sep="", collapse= ": ") 
  
  # if insert-genome junctions map to more than one base the "best" position is chosen 
  # based on the position with the most fragments, UMIs, reads
  # if two fragments have equal numbers of these both positions are kept
  
  if (length(unique(insert_reads$insert_end)) == 1) {
  insert_table$insert_base[i] <- unique(insert_reads$insert_end)
  } else {
    
    # print(cbind(unique(insert_reads$insert_end), frag_count, UMI_count, read_count))
    max_frags <- which(frag_count == max(as.numeric(unlist(frag_count)))) # highest fragment count
    max_UMIs  <- which(UMI_count  == max(as.numeric(unlist(UMI_count))))  # highest UMI count
    max_reads <- which(read_count == max(as.numeric(unlist(read_count)))) # highest read count 
    overlap   <- intersect(max_frags, max_UMIs)               # find overlap between frags & UMIs
    if(length(overlap) > 1) { overlap <- intersect(overlap, max_reads) } # if undecided try reads as a tiebreaker
    if(length(overlap) < 1) { overlap <- max_UMIs }                      # if no agreement between frags/UMIs/reads default to UMIs
    
    # ideally overlap is a single position but can be multiple positions if stats are identical
    insert_table$insert_base[i] <- toString(unique(insert_reads$insert_end)[overlap])
     
  }
}

# record a unique string for each insert combining sample, chr, strand and base
insert_table$insert_string  <- paste(insert_table$sample, insert_table$chr, insert_table$strand, insert_table$insert_base, sep=":")


# calculate clonality i.e. the relative abundance of each integration as a fraction of the total fragments/UMIs/reads
# calculate normalized clonality i.e. the clonality values normalized so that the most abundant value = 1
# clonality & normalized clonality allows quantitative comparison between samples with differing numbers of reads

insert_table$clon_fragments     <- insert_table$fragments/sum(insert_table$fragments)
insert_table$normclon_fragments <- insert_table$clon_fragments/max(insert_table$clon_fragments)

insert_table$clon_collUMIs      <- insert_table$collapsed_UMIs/sum(insert_table$collapsed_UMIs)
insert_table$normclon_collUMIs  <- insert_table$clon_collUMIs/max(insert_table$clon_collUMIs)

insert_table$clon_UMIs          <- insert_table$UMIs/sum(insert_table$UMIs)
insert_table$normclon_UMIs      <- insert_table$clon_UMIs/max(insert_table$clon_UMIs)

insert_table$clon_reads         <- insert_table$reads/sum(insert_table$reads)
insert_table$normclon_reads     <- insert_table$clon_reads/max(insert_table$clon_reads)



insert_table_long <- insert_table
# save the full table including UMIs and counts for each position
save(insert_table_long, file=paste("Tables/", sample_prefix, "_insert_table_long.RData", sep=""))

# save smaller table only including fields 1-18 that will be used in subsequent steps.
insert_table <- insert_table[,1:18]

save(insert_table, file=paste("Tables/", sample_prefix, "_insert_table.RData", sep=""))
write.table(insert_table, file=paste("Inserts/", sample_prefix, "_insert_table.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")


