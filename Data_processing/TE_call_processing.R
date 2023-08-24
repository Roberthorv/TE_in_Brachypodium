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

Full_TIP_TAP_Sample_table <- matrix(NA, nrow = dim(insertions)[1] + dim(deletions)[1], ncol = dim(coverage)[1])
colnames(Full_TIP_TAP_Sample_table) <- coverage[,1]
rownames(Full_TIP_TAP_Sample_table) <- c(insertions$ID, deletions$ID)

for (x in 1:dim(insertions)[1]) {
  Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==rownames(Full_TIP_TAP_Sample_table)[x], colnames(Full_TIP_TAP_Sample_table) %in% unlist(strsplit(insertions[insertions$ID==rownames(Full_TIP_TAP_Sample_table)[x],]$Presence_sample, split = ","))] <- 1
  Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==rownames(Full_TIP_TAP_Sample_table)[x], colnames(Full_TIP_TAP_Sample_table) %in% unlist(strsplit(insertions[insertions$ID==rownames(Full_TIP_TAP_Sample_table)[x],]$Absence_sample, split = ","))] <- 0
}
for (y in (dim(insertions)[1] +1):(dim(insertions)[1] + dim(deletions)[1])) {
  Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==rownames(Full_TIP_TAP_Sample_table)[y], colnames(Full_TIP_TAP_Sample_table) %in% unlist(strsplit(deletions[deletions$ID==rownames(Full_TIP_TAP_Sample_table)[y],]$Absence_sample, split = ","))] <- 1
  Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==rownames(Full_TIP_TAP_Sample_table)[y], colnames(Full_TIP_TAP_Sample_table) %in% unlist(strsplit(deletions[deletions$ID==rownames(Full_TIP_TAP_Sample_table)[y],]$Presence_sample, split = ","))] <- 0
}


# Make a TE info table
TE_info <- as.data.frame(matrix(NA, nrow = sum(dim(insertions)[1], dim(deletions)[1]), ncol = 6))
colnames(TE_info) <- c("TE_ID", "Scaffold", "Start_position", "End_position", "TE_name", "TE_family")
TE_info$TE_ID <- c(insertions$ID, deletions$ID)
TE_info$Scaffold <- c(insertions$Scaffold, deletions$Scaffold)
TE_info$Start_position <- c(insertions$Insertion_site_lower_estimate, deletions$Start_position)
TE_info$End_position <- c(insertions$Insertion_site_upper_estimate, deletions$End_position)
TE_info$TE_name <- c(insertions$TE_name, deletions$TE_name)
TE_info$TE_family <- c(insertions$TE_family, deletions$TE_family)
write.table(TE_info, file = "TE_info.txt", sep = "\t", header = TRUE)




