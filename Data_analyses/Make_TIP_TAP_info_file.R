## read in needed data
TE_info_consensus_family <- read.table("./TE_info_with_consensus_family.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Full_clean_TIP_TAP_Sample_table_342_samples <- read.table("./Full_clean_TIP_TAP_Sample_table_342_samples.txt", sep = "\t", header = TRUE)
Used_TIP_TAP_in_vcf <- Full_clean_TIP_TAP_Sample_table_342_samples[,c(1,2,4)]


TIP_TAP_in_vcf_info_bed <- Used_TIP_TAP_in_vcf[,c(1,2,2,3)]
TIP_TAP_in_vcf_info_bed[,2] <- TIP_TAP_in_vcf_info_bed[,2] - 1
TIP_TAP_in_vcf_info_bed[,5] <- NA
TIP_TAP_in_vcf_info_bed[,6] <- NA
TIP_TAP_in_vcf_info_bed[,7] <- NA


## add Class, Superfamily, family
for (i in 1:dim(TIP_TAP_in_vcf_info_bed)[1]) {
  my_ID <- unlist(strsplit(TIP_TAP_in_vcf_info_bed[i,4], split = "\\."))[1]
  my_info <- TE_info_consensus_family[TE_info_consensus_family[,1] == my_ID,]
  my_TE_sup_b <- unique(unlist(strsplit(as.character(my_info[6]), split = ",")))
  my_TE_sup <- my_TE_sup_b[my_TE_sup_b != "Unspecified"]
  if (length(my_TE_sup) < 1) {
    my_TE_sup <- NA
    my_TE_class <- NA
  } else{
    my_class <- c()
    for (j in 1:length(my_TE_sup)) {
      my_class <- c(my_class, unlist(strsplit(my_TE_sup[j], "/"))[1])
    }
    if (length(unique(my_class)) == 1) {
      if (unique(my_class) == "LTR") {
        my_TE_class <- "Retrotransposon"
      } else if (unique(my_class) == "DNA") {
        my_TE_class <- "DNA-transposon"
      } else if (unique(my_class) == "MITE") {
        my_TE_class <- "MITE(classIII)"
      }
    } else{my_TE_class <- NA}
    
  }
  
  TIP_TAP_in_vcf_info_bed[i,5] <- my_TE_class
  TIP_TAP_in_vcf_info_bed[i,6] <- paste(my_TE_sup, collapse = ",")
}

write.table(TIP_TAP_in_vcf_info_bed, file = "./TIP_TAP_in_vcf_info.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)




