# after runing TEPID
# read in the TEPID results
deletions <- as.data.frame(read.table("deletions_genotypes_without_ref.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(deletions) <- c("Scaffold", "Start_position", "End_position", "Strand", "TE_name", "Presence_sample", "Absence_sample")

insertions <- as.data.frame(read.table("insertions_genotypes_without_ref.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(insertions) <- c("Scaffold", "Insertion_site_lower_estimate", "Insertion_site_upper_estimate", "Donor_Scaffold", "Donor_Start_position", "Donor_End_position", "TE_name", "Presence_sample", "Absence_sample")

# read in the findings on the clade associations of each sample from Minadakis et al. 2023 
clade_table <- read.table("bdis_333_clade.csv", sep = ",", header = TRUE)

# read in a file with the name and coverage of all samples
coverage <- as.data.frame(read.table("Sample_coverage.txt",header = FALSE, sep="\t",stringsAsFactors=FALSE))
names(coverage) <- c("Sample", "Coverage")


# read in gff:
Bd_GFF <- as.data.frame(read.table("Bdistachyon_314_v3.0.fa.mod.EDTA.TEanno.gff3",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
Bd_GFF[,10:11] <- matrix(
  gsub("Classification=", "", 
       gsub("Name=", "", 
            unlist(strsplit(Bd_GFF[,9], split = ";"))[grepl(c("Name|Classification"),unlist(strsplit(Bd_GFF[,9], split = ";")))]
       )), ncol = 2, byrow = TRUE)

# Add TE family column:
insertions$TE_family <- NA
deletions$TE_family <- NA

for (k in 1:dim(insertions)[1]) {
  my_TE_names <- unlist(strsplit(insertions[k,]$TE_name, split = ","))
  my_TE_fam <- c()
  for (l in 1:length(my_TE_names)) {
    my_TE_fam <- c( my_TE_fam, unique(Bd_GFF[Bd_GFF[,10]==my_TE_names[l],11]))
  }
  if(length(my_TE_names)!=length(my_TE_fam)){break}
  else{insertions[k,]$TE_family <- paste(my_TE_fam, collapse = ",")}
}

for (k in 1:dim(deletions)[1]) {
  my_TE_names <- unlist(strsplit(deletions[k,]$TE_name, split = ","))
  my_TE_fam <- c()
  for (l in 1:length(my_TE_names)) {
    my_TE_fam <- c( my_TE_fam, unique(Bd_GFF[Bd_GFF[,10]==my_TE_names[l],11]))
  }
  if(length(my_TE_names)!=length(my_TE_fam)){break}
  else{deletions[k,]$TE_family <- paste(my_TE_fam, collapse = ",")}
}

# Add ID to the TIPs and TAPs:
insertions$ID <- paste0("TIP_",1:dim(insertions)[1])
deletions$ID <- paste0("TAP_",1:dim(deletions)[1])

## Make a Sample vs TIP/TAP table
# 1: TE is present; 0: TE is absent; NA no call

Full_TIP_TAP_Sample_table <- matrix(NA, nrow = dim(insertions)[1] + dim(deletions)[1], ncol = dim(coverage)[1] + 1)
colnames(Full_TIP_TAP_Sample_table) <- c(coverage[,1], "TE_ID")
rownames(Full_TIP_TAP_Sample_table) <- c(insertions$ID, deletions$ID)

for (x in 1:dim(insertions)[1]) {
  Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==rownames(Full_TIP_TAP_Sample_table)[x], colnames(Full_TIP_TAP_Sample_table) %in% unlist(strsplit(insertions[insertions$ID==rownames(Full_TIP_TAP_Sample_table)[x],]$Presence_sample, split = ","))] <- 1
  Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==rownames(Full_TIP_TAP_Sample_table)[x], colnames(Full_TIP_TAP_Sample_table) %in% unlist(strsplit(insertions[insertions$ID==rownames(Full_TIP_TAP_Sample_table)[x],]$Absence_sample, split = ","))] <- 0
}
for (y in (dim(insertions)[1] +1):(dim(insertions)[1] + dim(deletions)[1])) {
  Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==rownames(Full_TIP_TAP_Sample_table)[y], colnames(Full_TIP_TAP_Sample_table) %in% unlist(strsplit(deletions[deletions$ID==rownames(Full_TIP_TAP_Sample_table)[y],]$Absence_sample, split = ","))] <- 1
  Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==rownames(Full_TIP_TAP_Sample_table)[y], colnames(Full_TIP_TAP_Sample_table) %in% unlist(strsplit(deletions[deletions$ID==rownames(Full_TIP_TAP_Sample_table)[y],]$Presence_sample, split = ","))] <- 0
}
Full_TIP_TAP_Sample_table[,343] <- rownames(Full_TIP_TAP_Sample_table)

# Make a TE info table
TE_info <- as.data.frame(matrix(NA, nrow = sum(dim(insertions)[1], dim(deletions)[1]), ncol = 8))
colnames(TE_info) <- c("Scaffold", "Start_position", "End_position", "TE_ID", "TE_name", "TE_Superfamily", "Note_1", "Note_2")
TE_info$TE_ID <- c(insertions$ID, deletions$ID)
TE_info$Scaffold <- c(insertions$Scaffold, deletions$Scaffold)
TE_info$Start_position <- c(insertions$Insertion_site_lower_estimate, deletions$Start_position)
TE_info$End_position <- c(insertions$Insertion_site_upper_estimate, deletions$End_position)
TE_info$TE_name <- c(insertions$TE_name, deletions$TE_name)
TE_info$TE_family <- c(insertions$TE_family, deletions$TE_family)
TE_info$Note_1 <- NA
TE_info$Note_2 <- NA
write.table(TE_info, file = "TE_info.txt", sep = "\t", header = TRUE)

# merge tables
Full_TIP_TAP_Sample_table_with_ID_merge <- merge(x = TE_info, y = Full_TIP_TAP_Sample_table, by = "TE_ID")
Full_TIP_TAP_Sample_table_with_ID <- Full_TIP_TAP_Sample_table_with_ID_merge[,c(2,3,4,1,5:dim(Full_TIP_TAP_Sample_table_with_ID_merge)[2])]

# merge identical TE within 100 bp
TE_info_full_TIP_Sample_sorted <- Full_TIP_TAP_Sample_table_with_ID[grep("TIP",Full_TIP_TAP_Sample_table_with_ID[,4]),]

need_merge <- matrix(NA, ncol = 3, nrow = dim(TE_info_full_TIP_Sample_sorted)[1])
need_merge[,1] <- TE_info_full_TIP_Sample_sorted[,4]

need_merge[,2] <- unlist(sapply(1:dim(TE_info_full_TIP_Sample_sorted)[1], function(r){
  my_scaff <- TE_info_full_TIP_Sample_sorted[r,1]
  my_start_window <- TE_info_full_TIP_Sample_sorted[r,2] - 100
  my_end_window <- TE_info_full_TIP_Sample_sorted[r,3] + 100
  
  my_tab_for_merge <- TE_info_full_TIP_Sample_sorted[c(TE_info_full_TIP_Sample_sorted[,1] == my_scaff & TE_info_full_TIP_Sample_sorted[,2] >=  my_start_window & TE_info_full_TIP_Sample_sorted[,2] <= my_end_window) |
                                                           c(TE_info_full_TIP_Sample_sorted[,1] == my_scaff & TE_info_full_TIP_Sample_sorted[,3] >=  my_start_window & TE_info_full_TIP_Sample_sorted[,3] <= my_end_window),]
  
  if (length(my_tab_for_merge[my_tab_for_merge[,5] == TE_info_full_TIP_Sample_sorted[r,5],4]) == 1) {
    return("No")
  } else{
    return(paste(my_tab_for_merge[my_tab_for_merge[,5] == TE_info_full_TIP_Sample_sorted[r,5],4], collapse = ","))
  }
}))

need_merge[,3] <- unlist(sapply(1:dim(TE_info_full_TIP_Sample_sorted)[1], function(t){
  if (need_merge[t,2] == "No") {
    return(NA)
  } else{
    my_ID_to_merge <- unlist(strsplit(need_merge[t,2],split = ","))
    return(paste(unique(as.vector(unlist(sapply(my_ID_to_merge, function(s){
      my_cand <- grep(pattern = s, need_merge[,2])
      true_line_num <- c()
      for (u in 1:length(my_cand)) {
        if (length(which(unlist(strsplit(need_merge[my_cand[u],2], split = ",")) == s))==1) {
          true_line_num <- c(true_line_num, my_cand[u])
        } else next
      }
      return(true_line_num)
    })))), collapse = ","))
  }
}))

lines_to_merge <- unique(need_merge[!is.na(need_merge[,3]),3])

merged_TE_calls <- matrix(NA, nrow = length(lines_to_merge), ncol = (1 + dim(TE_info_full_TIP_Sample_sorted)[2]))
colnames(merged_TE_calls) <- c(colnames(TE_info_full_TIP_Sample_sorted[1:4]), "merged_TE_ID",colnames(TE_info_full_TIP_Sample_sorted)[5:350])

for (l in 1:length(lines_to_merge)) {
  my_lines_full <- TE_info_full_TIP_Sample_sorted[as.numeric(unlist(strsplit(lines_to_merge[l], split = ","))),]
  merged_TE_calls[l,1:3] <- unlist(my_lines_full[1,1:3])
  merged_TE_calls[l,4] <- paste0(my_lines_full[1,4], ".1")
  merged_TE_calls[l,5] <- paste(my_lines_full[,4], collapse = ",")
  merged_TE_calls[l,6:9] <- unlist(my_lines_full[1,5:8])
  my_lines_full_num <- as.data.frame(my_lines_full[,9:350])
  for  (o in 1:dim(my_lines_full_num)[2]){
    my_lines_full_num[,o] <- as.numeric(my_lines_full_num[,o])
  }
  counts <- colSums(my_lines_full_num)
  counts[counts>1 & !is.na(counts)] <- 1
  merged_TE_calls[l,10:351] <- counts
}

TE_info_full_TAP_Sample_sorted <- Full_TIP_TAP_Sample_table_with_ID[grep("TAP",Full_TIP_TAP_Sample_table_with_ID[,4]),]

TAP_no_need_for_merging_TE_calls <- matrix(NA, nrow = dim(TE_info_full_TAP_Sample_sorted)[1], ncol = (1 + dim(TE_info_full_TAP_Sample_sorted)[2]))
colnames(TAP_no_need_for_merging_TE_calls) <- c(colnames(TE_info_full_TAP_Sample_sorted[1:4]), "merged_TE_ID",colnames(TE_info_full_TAP_Sample_sorted)[5:350])
TAP_no_need_for_merging_TE_calls[,1:4] <- unlist(TE_info_full_TAP_Sample_sorted[,1:4])
TAP_no_need_for_merging_TE_calls[,6:351] <- unlist(TE_info_full_TAP_Sample_sorted[,5:350])

no_need_for_merging_TE_calls <- matrix(NA, nrow = dim(TE_info_full_TIP_Sample_sorted[need_merge[,2] == "No",])[1], ncol = (1 + dim(TE_info_full_TIP_Sample_sorted)[2]))
colnames(no_need_for_merging_TE_calls) <- c(colnames(TE_info_full_TIP_Sample_sorted[1:4]), "merged_TE_ID",colnames(TE_info_full_TIP_Sample_sorted)[5:350])
no_need_for_merging_TE_calls[,1:4] <- unlist(TE_info_full_TIP_Sample_sorted[need_merge[,2] == "No",1:4])
no_need_for_merging_TE_calls[,6:351] <- unlist(TE_info_full_TIP_Sample_sorted[need_merge[,2] == "No",5:350])

Full_merged_TE_calls <- rbind(TAP_no_need_for_merging_TE_calls, no_need_for_merging_TE_calls, merged_TE_calls)
Full_merged_TE_calls_df <- as.data.frame(Full_merged_TE_calls)
write.table(Full_merged_TE_calls, file = "Full_merged_TE_calls.txt", sep = "\t", quote = FALSE, row.names = FALSE)











