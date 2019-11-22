# if a sample is absent then it had no mappable inserts (ideal for control samples)
insertfilesfull <- dir(path = "Tables/", pattern = "insert_table_long.RData")
insertfilesfull <- paste("Tables/",insertfilesfull,sep="")

load(file = insertfilesfull[1])
print(insertfilesfull[1])
all_inserts_full <- insert_table_long

# rbind each table to all_inserts, this loop gets quite slow when > 500,000 lines due to use of rbind
# if needed can break the loop up into multiple sets of insert files and do rbind of multiple all_insert files at the end

for (insertfilefull in 2:length(insertfilesfull)){
  print(insertfilesfull[insertfilefull])
  load(file = insertfilesfull[insertfilefull])
  all_inserts_full <- rbind(all_inserts_full, insert_table_long)
}

# insert_base is still recorded as a string that occasionally has multiple values
# grep(",", all_inserts$insert_base)
# to identify duplicates record insert base as a numeric with an average value of multiple entries

# an alternate strategy (not implemented here) would be to insert_table_long to record multiple lines that are identical 
# except for the bases, these can later be regrouped

average_base <- numeric()
for (i in 1:nrow(all_inserts_full)) {
  average_base[i] <- mean(as.numeric(unlist(strsplit(all_inserts_full$insert_base[i],", "))))
}
all_inserts_full <- cbind(all_inserts_full, average_base, stringsAsFactors = FALSE)

# add a column concatenating chromosome and strand for batch processing of the next step
chr_strand <- paste(all_inserts_full$chr, all_inserts_full$strand,sep=":")
all_inserts_full <- cbind(all_inserts_full, chr_strand, stringsAsFactors = FALSE)

# order by this column
all_inserts_full <- all_inserts_full[order(all_inserts_full$chr_strand, all_inserts_full$average_base), ]

# add insert_id and empty columns to record duplicates and the categories of duplicates
insert_id        <- 1:nrow(all_inserts_full)
duplicates       <- rep("",nrow(all_inserts_full))
duplicate_coords <- rep("",nrow(all_inserts_full))
control          <- rep("",nrow(all_inserts_full))
replicate_groups <- rep("",nrow(all_inserts_full))
expected         <- rep("",nrow(all_inserts_full))
rescue           <- rep("",nrow(all_inserts_full))

row.names(all_inserts_full) <- insert_id
all_inserts_full<- cbind(all_inserts_full, duplicates, duplicate_coords, control, replicate_groups, expected, rescue, stringsAsFactors = FALSE)



