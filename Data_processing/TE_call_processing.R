# after runing TEPID
# read in the TEPID results
deletions <- as.data.frame(read.table("deletions_genotypes_without_ref.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(deletions) <- c("Scaffold", "Start_position", "End_position", "Strand", "TE_name", "Presence_sample", "Absence_sample")

insertions <- as.data.frame(read.table("insertions_genotypes_without_ref.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(insertions) <- c("Scaffold", "Insertion_site_lower_estimate", "Insertion_site_upper_estimate", "Donor_Scaffold", "Donor_Start_position", "Donor_End_position", "TE_name", "Presence_sample", "Absence_sample")

# read in the findings on the clade associations of each sample from Minadakis et al. 2023 
clade_table <- read.table("bdis_333_clade.csv", sep = ",", header = TRUE)



