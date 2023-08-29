#!/usr/bin/env Rscript

args <-   commandArgs(trailingOnly=TRUE)
library(vcfR, lib.loc = "./R_packages/")
library(memuse, lib.loc = "./R_packages/")

print(args)

vcf <- read.vcfR("./small_mock_vcf_for_TEs.vcf", verbose = FALSE )
Full_clean_TIP_TAP_Sample_table_342_samples <- read.table("./Full_clean_TIP_TAP_Sample_table_342_samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

for (i in args[1]:args[2]) {
  my_chrom <- Full_clean_TIP_TAP_Sample_table_342_samples[i,1]
  my_ID <- Full_clean_TIP_TAP_Sample_table_342_samples[i,4]
  my_qual <- 1000
  my_filter <- "PASS"
  my_AC <- sum(Full_clean_TIP_TAP_Sample_table_342_samples[i, 11:dim(Full_clean_TIP_TAP_Sample_table_342_samples)[2]], na.rm = TRUE)
  my_AN <- sum(!is.na(Full_clean_TIP_TAP_Sample_table_342_samples[i, 11:dim(Full_clean_TIP_TAP_Sample_table_342_samples)[2]]))
  my_AF <- format(my_AC/my_AN, scientific = TRUE)
  my_info <- paste(c("AC=", my_AC, ";AF=", my_AF, ";AN=", my_AN), collapse = "")
  if (grepl(",", Full_clean_TIP_TAP_Sample_table_342_samples[i,6])) {
    multiallelic <- TRUE
    my_info <- paste0(my_info, ";MULTIALLELIC")
  } else{multiallelic <- FALSE}
  if (!multiallelic) {
    if (grepl("TIP", Full_clean_TIP_TAP_Sample_table_342_samples[i,4])) {
      my_pos <- ceiling((Full_clean_TIP_TAP_Sample_table_342_samples[i,2] + Full_clean_TIP_TAP_Sample_table_342_samples[i,3])/2)
      my_ref <- "A"
      my_alt <- "T"
    } else{
      my_pos <- Full_clean_TIP_TAP_Sample_table_342_samples[i,2]
      my_ref <- "T"
      my_alt <- "A"
    }
  } else{
    my_pos <- ceiling((Full_clean_TIP_TAP_Sample_table_342_samples[i,2] + Full_clean_TIP_TAP_Sample_table_342_samples[i,3])/2)
    my_ref <- "A"
    my_n_alt <- length(unlist(strsplit(Full_clean_TIP_TAP_Sample_table_342_samples[i,4], split = ",")))
    my_alt <- "T"
    for (h in 2:my_n_alt) {
      my_alt <- paste(my_alt, paste0(c("T", rep("A", h-1)), collapse = ""), sep = ",")
    }
  }
  
  vcf@fix[1 + (i - as.numeric(args[1])),1] <- my_chrom
  vcf@fix[1 + (i - as.numeric(args[1])),2] <- my_pos
  vcf@fix[1 + (i - as.numeric(args[1])),3] <- my_ID
  vcf@fix[1 + (i - as.numeric(args[1])),4] <- my_ref
  vcf@fix[1 + (i - as.numeric(args[1])),5] <- my_alt
  vcf@fix[1 + (i - as.numeric(args[1])),6] <- my_qual
  vcf@fix[1 + (i - as.numeric(args[1])),7] <- my_filter
  vcf@fix[1 + (i - as.numeric(args[1])),8] <- my_info
  vcf@gt[1 + (i - as.numeric(args[1])),1] <- "GT"
  
  for (x in 11:dim(Full_clean_TIP_TAP_Sample_table_342_samples)[2]) {
    if (is.na(Full_clean_TIP_TAP_Sample_table_342_samples[i,x])) {
      vcf@gt[1 + (i - as.numeric(args[1])),1 + (x -10)] <- "./."
    } else{
      if (Full_clean_TIP_TAP_Sample_table_342_samples[i,x] == 1) {
        vcf@gt[1 + (i - as.numeric(args[1])),1 + (x -10)] <- "1/1"
      } else{vcf@gt[1 + (i - as.numeric(args[1])),1 + (x -10)] <- "0/0"}}
  }
}

write.vcf(vcf, paste0("Full_clean_TE_call_TEPID_", args[1], "_", args[2], ".vcf.gz"))
