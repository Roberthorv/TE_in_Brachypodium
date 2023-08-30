## Make an age adjusted SFS for each clade
## use a freqeuncy adjustment at the begining to account for frequency call biases 

library(vcfR)

## get Retrotransposon 
TIP_TAP_in_vcf_info_bed <- read.table("./TIP_TAP_in_vcf_info.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Retrotransposon_ID <- TIP_TAP_in_vcf_info_bed[TIP_TAP_in_vcf_info_bed[,5] == "Retrotransposon" & !is.na(TIP_TAP_in_vcf_info_bed[,5]), 4]
