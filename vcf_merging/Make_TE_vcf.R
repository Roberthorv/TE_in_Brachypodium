#!/usr/bin/env Rscript

args <-   commandArgs(trailingOnly=TRUE)
library(vcfR, lib.loc = "./R_packages/")
library(memuse, lib.loc = "./R_packages/")

vcf <- read.vcfR("./mock_vcf_for_TEs.vcf", verbose = FALSE )
TE_info <- read.table("./TE_info.txt", header = TRUE, stringsAsFactors = FALSE)
Full_TIP_TAP_Sample_table <- read.table("./Full_clean_TIP_TAP_Sample_table_342_samples.txt", header = TRUE, stringsAsFactors = FALSE)

for (i in args[1]:args[2]) {
  my_chrom <- TE_info[i,]$Scaffold
  my_ID <- TE_info[i,]$TE_ID
  my_qual <- 1000
  my_filter <- "PASS"
  my_AC <- sum(Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==TE_info[i,]$TE_ID, ], na.rm = TRUE)
  my_AN <- sum(!is.na(Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==TE_info[i,]$TE_ID, ]))
  my_AF <- format(my_AC/my_AN, scientific = TRUE)
  my_info <- paste(c("AC=", my_AC, ";AF=", my_AF, ";AN=", my_AN), collapse = "")
  if (grepl(",", TE_info[i,]$TE_name)) {
    multiallelic <- TRUE
    my_info <- paste0(my_info, ";MULTIALLELIC")
  } else{multiallelic <- FALSE}
  if (!multiallelic) {
    if (grepl("TIP", TE_info[i,]$TE_ID)) {
      my_pos <- ceiling((TE_info[i,]$Start_position + TE_info[i,]$End_position)/2)
      my_ref <- "A"
      my_alt <- "T"
    }
    else{
      my_pos <- TE_info[i,]$Start_position
      my_ref <- "T"
      my_alt <- "A"
    }
  }
  else{
    my_pos <- ceiling((TE_info[i,]$Start_position + TE_info[i,]$End_position)/2)
    my_ref <- "A"
    my_n_alt <- length(unlist(strsplit(TE_info[i,]$TE_name, split = ",")))
    my_alt <- "T"
    for (h in 2:my_n_alt) {
      my_alt <- paste(my_alt, paste0(c("T", rep("A", h-1)), collapse = ""), sep = ",")
    }
  }
  
  vcf@fix[i,1] <- my_chrom
  vcf@fix[i,2] <- my_pos
  vcf@fix[i,3] <- my_ID
  vcf@fix[i,4] <- my_ref
  vcf@fix[i,5] <- my_alt
  vcf@fix[i,6] <- my_qual
  vcf@fix[i,7] <- my_filter
  vcf@fix[i,8] <- my_info
  vcf@gt[i,1] <- "GT"
  
  for (x in 1:dim(Full_TIP_TAP_Sample_table)[2]) {
    if (is.na(Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==TE_info[i,]$TE_ID, x])) {
      vcf@gt[i,1 + x] <- "./."
    } else{
      if (Full_TIP_TAP_Sample_table[rownames(Full_TIP_TAP_Sample_table)==TE_info[i,]$TE_ID, x] == 1) {
        vcf@gt[i,1 + x] <- "1/1"
      } else{vcf@gt[i,1 + x] <- "0/0"}}
  }
}

write.vcf(vcf, paste0("Full_TE_call_TEPID_", args[1], "_", args[2], ".vcf"))


