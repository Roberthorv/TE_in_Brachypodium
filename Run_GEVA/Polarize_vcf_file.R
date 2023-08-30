#!/usr/bin/env Rscript

args <-   commandArgs(trailingOnly=TRUE)
library(vcfR, lib.loc = "./R_packages/")
library(memuse, lib.loc = "./R_packages/")

vcf <- read.vcfR(paste0("./To_fixed_SNPs_TE_part_", args[1], ".vcf"), verbose = FALSE)

for (i in 1:dim(vcf@fix)[1]) {
  print(i/dim(vcf@fix)[1])
  my_anc <- vcf@fix[i,5]
  my_der <- vcf@fix[i,4]
  my_info <- vcf@fix[i,8]
  my_info_to_replace <- paste(unlist(strsplit(my_info, split = ";"))[1:3], collapse = ";")
  my_AC <- paste("AC=", (as.numeric(gsub(pattern = "AN=", replacement = "",unlist(strsplit(my_info_to_replace, split = ";"))[3])) - as.numeric(gsub(pattern = "AC=", replacement = "",unlist(strsplit(my_info_to_replace, split = ";"))[1]))), sep = "")
  my_AF <- paste("AF=", (1- as.numeric(gsub(pattern = "AF=", replacement = "",unlist(strsplit(my_info_to_replace, split = ";"))[2])) ), sep = "")
  my_info_replace <- paste(my_AC, my_AF, unlist(strsplit(my_info_to_replace, split = ";"))[3], sep = ";")
  
  vcf@fix[i,4] <- my_anc
  vcf@fix[i,5] <- my_der
  vcf@fix[i,8] <- gsub(pattern = my_info_to_replace , replacement = my_info_replace, vcf@fix[i,8])
  vcf@gt[i,] <- gsub(pattern = "X/X" , replacement = "0/0", gsub(pattern = "0/0" , replacement = "1/1", gsub(pattern = "1/1" , replacement = "X/X", vcf@gt[i,])))

}

write.vcf(vcf, paste0("./Pol_SNPs_TE_part_", args[1], ".vcf.gz"))

