# if a sample is absent then it had no mappable inserts (ideal for control samples)
insertfiles <- dir(path = "Tables/", pattern = "insert_table.RData")
insertfiles <- paste("Tables/",insertfiles,sep="")

load(file = insertfiles[1])
print(insertfiles[1])
all_inserts <- insert_table

# rbind each table to all_inserts, this loop gets quite slow when > 500,000 lines due to use of rbind
# if needed can break the loop up into multiple sets of insert files and do rbind of multiple all_insert files at the end

for (insertfile in 2:length(insertfiles)){
  print(insertfiles[insertfile])
  load(file = insertfiles[insertfile])
  all_inserts <- rbind(all_inserts, insert_table) 
}

# insert_base is still recorded as a string that occasionally has multiple values
# grep(",", all_inserts$insert_base)
# to identify duplicates record insert base as a numeric with an average value of multiple entries

# an alternate strategy (not implemented here) would be to insert_table_long to record multiple lines that are identical 
# except for the bases, these can later be regrouped

average_base <- numeric()
for (i in 1:nrow(all_inserts)) {
  average_base[i] <- mean(as.numeric(unlist(strsplit(all_inserts$insert_base[i],", "))))
}
all_inserts <- cbind(all_inserts, average_base, stringsAsFactors = FALSE)

# add a column concatenating chromosome and strand for batch processing of the next step
chr_strand <- paste(all_inserts$chr, all_inserts$strand,sep=":")
all_inserts <- cbind(all_inserts, chr_strand, stringsAsFactors = FALSE)

# order by this column
all_inserts <- all_inserts[order(all_inserts$chr_strand, all_inserts$average_base), ]

# add insert_id and empty columns to record duplicates and the categories of duplicates
insert_id        <- 1:nrow(all_inserts)
duplicates       <- rep("",nrow(all_inserts))
duplicate_coords <- rep("",nrow(all_inserts))
control          <- rep("",nrow(all_inserts))
replicate_groups <- rep("",nrow(all_inserts))
expected         <- rep("",nrow(all_inserts))
rescue           <- rep("",nrow(all_inserts))

row.names(all_inserts) <- insert_id
all_inserts<- cbind(insert_id, all_inserts, duplicates, duplicate_coords, control, replicate_groups, expected, rescue, stringsAsFactors = FALSE)
