## read in TIP TAP data:
Full_clean_TIP_TAP_Sample_table_342_samples <- read.table("Full_clean_TIP_TAP_Sample_table_342_samples.txt", sep = "\t", header = TRUE)

## get the 326 samples:
## some samples were removed because there were outliers in the PCA (Bd1-1 Bd21-3 BdTR8i BdTR3C Bd18-1 Bd21)
Full_clean_TIP_TAP_Sample_table_326_samples_all <- Full_clean_TIP_TAP_Sample_table_342_samples[,!names(Full_clean_TIP_TAP_Sample_table_342_samples) %in% c("ABR8", "Bd21Control_1", "Bd21Control_2", "Bd21Control_3", "Bd21_Skalska", "Cef2", "D26", "D73", "Pob1", "SRR1800522", "Bd21_Stritt", "Bd1.1", "Bd21.3", "BdTR8i", "BdTR3C", "Bd18.1")]

row_to_exclude <- sapply(1:dim(Full_clean_TIP_TAP_Sample_table_326_samples_all)[1], function(x){
  n_allele_1 <- sum(Full_clean_TIP_TAP_Sample_table_326_samples_all[x, 11:336] == 0, na.rm = TRUE)
  n_allele_2 <- sum(Full_clean_TIP_TAP_Sample_table_326_samples_all[x, 11:336] == 1, na.rm = TRUE)
  if (n_allele_2 * n_allele_1 == 0) {
    return(x)
  }
})

Full_clean_TIP_TAP_Sample_table_326_samples <- Full_clean_TIP_TAP_Sample_table_326_samples_all[-unlist(row_to_exclude),]
write.table(Full_clean_TIP_TAP_Sample_table_326_samples, file = "Full_clean_TIP_TAP_Sample_table_326_samples.txt", quote = FALSE, sep = "\t")

## read in TIP TAP data 326 samples:
Full_clean_TIP_TAP_Sample_table_326_samples <- read.table("Full_clean_TIP_TAP_Sample_table_326_samples.txt", sep = "\t", header = TRUE)

## Total number of polymorphic TE insertions 
dim(Full_clean_TIP_TAP_Sample_table_326_samples)[1]

## Total number of retrotransposons and DNA transposons 
TIP_TAP_in_vcf_info_bed <- read.table("TIP_TAP_in_vcf_info.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

## get Retrotransposon 
Retrotransposon_ID <- TIP_TAP_in_vcf_info_bed[TIP_TAP_in_vcf_info_bed[,5] == "Retrotransposon" & !is.na(TIP_TAP_in_vcf_info_bed[,5]), 4]
dim(Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples$TE_ID %in% Retrotransposon_ID,])[1]

## get DNAtransposon
DNAtransposon_ID <- TIP_TAP_in_vcf_info_bed[TIP_TAP_in_vcf_info_bed[,5] == "DNA-transposon" & !is.na(TIP_TAP_in_vcf_info_bed[,5]), 4]
dim(Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples$TE_ID %in% DNAtransposon_ID,])[1]

## run a PCA:
library(devtools)
library(ggbiplot)
library(vcfR)
library(ggpubr)

## run PCA on SNPs

## read in data:
four_fold_degenerate_SNP_Bdis326_vcf <- read.vcfR("four_fold_degenerate_SNP_Bdis326.vcf.gz", verbose = FALSE )

## extract SNPs
four_fold_degenerate_SNP_Bdis326_extract_gt <- extract.gt(four_fold_degenerate_SNP_Bdis326_vcf, element = "GT", IDtoRowNames  = F, as.numeric = T, convertNA = T, return.alleles = F)
four_fold_degenerate_SNP_Bdis326_extract_gt_df <- data.frame(t(four_fold_degenerate_SNP_Bdis326_extract_gt[complete.cases(four_fold_degenerate_SNP_Bdis326_extract_gt),]))

df_pca_data_four_fold_degenerate_SNP_Bdis326 <- four_fold_degenerate_SNP_Bdis326_extract_gt_df[,apply(four_fold_degenerate_SNP_Bdis326_extract_gt_df, 2, var, na.rm=TRUE) != 0]

pca_four_fold_degenerate_SNP_Bdis326.pca <- prcomp(df_pca_data_four_fold_degenerate_SNP_Bdis326, center = T, scale. = T)





