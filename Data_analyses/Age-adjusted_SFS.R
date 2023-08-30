## Make an age adjusted SFS for each clade
## use a freqeuncy adjustment at the begining to account for frequency call biases 

library(vcfR)

## get Retrotransposon 
TIP_TAP_in_vcf_info_bed <- read.table("./TIP_TAP_in_vcf_info.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Retrotransposon_ID <- TIP_TAP_in_vcf_info_bed[TIP_TAP_in_vcf_info_bed[,5] == "Retrotransposon" & !is.na(TIP_TAP_in_vcf_info_bed[,5]), 4]

## read in the results from geva

## A Italia
{
  synonymous_A_Italia.vcf <- read.vcfR("./VCF_files/four_fold_degenerate_SNP_Bdis326_A_Italia.vcf", verbose = FALSE )
  
  Bd1_A_Italia_Clock_J_age_estimate <- read.table("./GEVA/A_Italia/Bd1_A_Italia_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd2_A_Italia_Clock_J_age_estimate <- read.table("./GEVA/A_Italia/Bd2_A_Italia_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd3_A_Italia_Clock_J_age_estimate <- read.table("./GEVA/A_Italia/Bd3_A_Italia_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd4_A_Italia_Clock_J_age_estimate <- read.table("./GEVA/A_Italia/Bd4_A_Italia_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd5_A_Italia_Clock_J_age_estimate <- read.table("./GEVA/A_Italia/Bd5_A_Italia_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd1_A_Italia.marker <- read.table("./GEVA/A_Italia/Bd1_A_Italia.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd2_A_Italia.marker <- read.table("./GEVA/A_Italia/Bd2_A_Italia.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd3_A_Italia.marker <- read.table("./GEVA/A_Italia/Bd3_A_Italia.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd4_A_Italia.marker <- read.table("./GEVA/A_Italia/Bd4_A_Italia.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd5_A_Italia.marker <- read.table("./GEVA/A_Italia/Bd5_A_Italia.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  not_derived_TE_Bd1_A_Italia <- read.table("./GEVA/A_Italia/not_derived_TE_Bd1_A_Italia.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd2_A_Italia <- read.table("./GEVA/A_Italia/not_derived_TE_Bd2_A_Italia.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd3_A_Italia <- read.table("./GEVA/A_Italia/not_derived_TE_Bd3_A_Italia.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd4_A_Italia <- read.table("./GEVA/A_Italia/not_derived_TE_Bd4_A_Italia.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd5_A_Italia <- read.table("./GEVA/A_Italia/not_derived_TE_Bd5_A_Italia.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  not_derived_TE_All_A_Italia <- rbind(not_derived_TE_Bd1_A_Italia, not_derived_TE_Bd2_A_Italia, not_derived_TE_Bd3_A_Italia, not_derived_TE_Bd4_A_Italia, not_derived_TE_Bd5_A_Italia)  
  
  Bd1_A_Italia.marker$Frequencies <- Bd1_A_Italia.marker$AlleleCount1/(Bd1_A_Italia.marker$AlleleCount1 + Bd1_A_Italia.marker$AlleleCount0)
  Bd1_A_Italia.marker$Sigelton <- unlist(sapply(1:dim(Bd1_A_Italia.marker)[1], function(x){
    if (Bd1_A_Italia.marker$AlleleCount1[x] == 2 | Bd1_A_Italia.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd1_A_Italia_merge <- merge(Bd1_A_Italia.marker, Bd1_A_Italia_Clock_J_age_estimate, by = "MarkerID")
  
  Bd2_A_Italia.marker$Frequencies <- Bd2_A_Italia.marker$AlleleCount1/(Bd2_A_Italia.marker$AlleleCount1 + Bd2_A_Italia.marker$AlleleCount0)
  Bd2_A_Italia.marker$Sigelton <- unlist(sapply(1:dim(Bd2_A_Italia.marker)[1], function(x){
    if (Bd2_A_Italia.marker$AlleleCount1[x] == 2 | Bd2_A_Italia.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd2_A_Italia_merge <- merge(Bd2_A_Italia.marker, Bd2_A_Italia_Clock_J_age_estimate, by = "MarkerID")
  
  Bd3_A_Italia.marker$Frequencies <- Bd3_A_Italia.marker$AlleleCount1/(Bd3_A_Italia.marker$AlleleCount1 + Bd3_A_Italia.marker$AlleleCount0)
  Bd3_A_Italia.marker$Sigelton <- unlist(sapply(1:dim(Bd3_A_Italia.marker)[1], function(x){
    if (Bd3_A_Italia.marker$AlleleCount1[x] == 2 | Bd3_A_Italia.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd3_A_Italia_merge <- merge(Bd3_A_Italia.marker, Bd3_A_Italia_Clock_J_age_estimate, by = "MarkerID")
  
  Bd4_A_Italia.marker$Frequencies <- Bd4_A_Italia.marker$AlleleCount1/(Bd4_A_Italia.marker$AlleleCount1 + Bd4_A_Italia.marker$AlleleCount0)
  Bd4_A_Italia.marker$Sigelton <- unlist(sapply(1:dim(Bd4_A_Italia.marker)[1], function(x){
    if (Bd4_A_Italia.marker$AlleleCount1[x] == 2 | Bd4_A_Italia.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd4_A_Italia_merge <- merge(Bd4_A_Italia.marker, Bd4_A_Italia_Clock_J_age_estimate, by = "MarkerID")
  
  Bd5_A_Italia.marker$Frequencies <- Bd5_A_Italia.marker$AlleleCount1/(Bd5_A_Italia.marker$AlleleCount1 + Bd5_A_Italia.marker$AlleleCount0)
  Bd5_A_Italia.marker$Sigelton <- unlist(sapply(1:dim(Bd5_A_Italia.marker)[1], function(x){
    if (Bd5_A_Italia.marker$AlleleCount1[x] == 2 | Bd5_A_Italia.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd5_A_Italia_merge <- merge(Bd5_A_Italia.marker, Bd5_A_Italia_Clock_J_age_estimate, by = "MarkerID")
  
  All_A_Italia_merge <- rbind(Bd1_A_Italia_merge, Bd2_A_Italia_merge, Bd3_A_Italia_merge, Bd4_A_Italia_merge, Bd5_A_Italia_merge)
  
  synonymous_Bd1_A_Italia_position <- as.numeric(synonymous_A_Italia.vcf@fix[synonymous_A_Italia.vcf@fix[,1]=="Bd1",2])
  synonymous_Bd1_A_Italia_age <- Bd1_A_Italia_merge[Bd1_A_Italia_merge[,3]%in%synonymous_Bd1_A_Italia_position, c(1:3, 15, 20)]
  
  synonymous_Bd2_A_Italia_position <- as.numeric(synonymous_A_Italia.vcf@fix[synonymous_A_Italia.vcf@fix[,1]=="Bd2",2])
  synonymous_Bd2_A_Italia_age <- Bd2_A_Italia_merge[Bd2_A_Italia_merge[,3]%in%synonymous_Bd2_A_Italia_position, c(1:3, 15, 20)]
  
  synonymous_Bd3_A_Italia_position <- as.numeric(synonymous_A_Italia.vcf@fix[synonymous_A_Italia.vcf@fix[,1]=="Bd3",2])
  synonymous_Bd3_A_Italia_age <- Bd3_A_Italia_merge[Bd3_A_Italia_merge[,3]%in%synonymous_Bd3_A_Italia_position, c(1:3, 15, 20)]
  
  synonymous_Bd4_A_Italia_position <- as.numeric(synonymous_A_Italia.vcf@fix[synonymous_A_Italia.vcf@fix[,1]=="Bd4",2])
  synonymous_Bd4_A_Italia_age <- Bd4_A_Italia_merge[Bd4_A_Italia_merge[,3]%in%synonymous_Bd4_A_Italia_position, c(1:3, 15, 20)]
  
  synonymous_Bd5_A_Italia_position <- as.numeric(synonymous_A_Italia.vcf@fix[synonymous_A_Italia.vcf@fix[,1]=="Bd5",2])
  synonymous_Bd5_A_Italia_age <- Bd5_A_Italia_merge[Bd5_A_Italia_merge[,3]%in%synonymous_Bd5_A_Italia_position, c(1:3, 15, 20)]
  
  synonymous_All_A_Italia_age <- rbind(synonymous_Bd1_A_Italia_age, synonymous_Bd2_A_Italia_age, synonymous_Bd3_A_Italia_age, synonymous_Bd4_A_Italia_age, synonymous_Bd5_A_Italia_age)
  
}
## A East
{
  synonymous_A_East.vcf <- read.vcfR("./VCF_files/four_fold_degenerate_SNP_Bdis326_A_East.vcf", verbose = FALSE )
  
  Bd1_A_East_Clock_J_age_estimate <- read.table("./GEVA/A_East/Bd1_A_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd2_A_East_Clock_J_age_estimate <- read.table("./GEVA/A_East/Bd2_A_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd3_A_East_Clock_J_age_estimate <- read.table("./GEVA/A_East/Bd3_A_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd4_A_East_Clock_J_age_estimate <- read.table("./GEVA/A_East/Bd4_A_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd5_A_East_Clock_J_age_estimate <- read.table("./GEVA/A_East/Bd5_A_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd1_A_East.marker <- read.table("./GEVA/A_East/Bd1_A_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd2_A_East.marker <- read.table("./GEVA/A_East/Bd2_A_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd3_A_East.marker <- read.table("./GEVA/A_East/Bd3_A_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd4_A_East.marker <- read.table("./GEVA/A_East/Bd4_A_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd5_A_East.marker <- read.table("./GEVA/A_East/Bd5_A_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  not_derived_TE_Bd1_A_East <- read.table("./GEVA/A_East/not_derived_TE_Bd1_A_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd2_A_East <- read.table("./GEVA/A_East/not_derived_TE_Bd2_A_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd3_A_East <- read.table("./GEVA/A_East/not_derived_TE_Bd3_A_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd4_A_East <- read.table("./GEVA/A_East/not_derived_TE_Bd4_A_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd5_A_East <- read.table("./GEVA/A_East/not_derived_TE_Bd5_A_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  not_derived_TE_All_A_East <- rbind(not_derived_TE_Bd1_A_East, not_derived_TE_Bd2_A_East, not_derived_TE_Bd3_A_East, not_derived_TE_Bd4_A_East, not_derived_TE_Bd5_A_East)  
  
  Bd1_A_East.marker$Frequencies <- Bd1_A_East.marker$AlleleCount1/(Bd1_A_East.marker$AlleleCount1 + Bd1_A_East.marker$AlleleCount0)
  Bd1_A_East.marker$Sigelton <- unlist(sapply(1:dim(Bd1_A_East.marker)[1], function(x){
    if (Bd1_A_East.marker$AlleleCount1[x] == 2 | Bd1_A_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd1_A_East_merge <- merge(Bd1_A_East.marker, Bd1_A_East_Clock_J_age_estimate, by = "MarkerID")
  
  Bd2_A_East.marker$Frequencies <- Bd2_A_East.marker$AlleleCount1/(Bd2_A_East.marker$AlleleCount1 + Bd2_A_East.marker$AlleleCount0)
  Bd2_A_East.marker$Sigelton <- unlist(sapply(1:dim(Bd2_A_East.marker)[1], function(x){
    if (Bd2_A_East.marker$AlleleCount1[x] == 2 | Bd2_A_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd2_A_East_merge <- merge(Bd2_A_East.marker, Bd2_A_East_Clock_J_age_estimate, by = "MarkerID")
  
  Bd3_A_East.marker$Frequencies <- Bd3_A_East.marker$AlleleCount1/(Bd3_A_East.marker$AlleleCount1 + Bd3_A_East.marker$AlleleCount0)
  Bd3_A_East.marker$Sigelton <- unlist(sapply(1:dim(Bd3_A_East.marker)[1], function(x){
    if (Bd3_A_East.marker$AlleleCount1[x] == 2 | Bd3_A_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd3_A_East_merge <- merge(Bd3_A_East.marker, Bd3_A_East_Clock_J_age_estimate, by = "MarkerID")
  
  Bd4_A_East.marker$Frequencies <- Bd4_A_East.marker$AlleleCount1/(Bd4_A_East.marker$AlleleCount1 + Bd4_A_East.marker$AlleleCount0)
  Bd4_A_East.marker$Sigelton <- unlist(sapply(1:dim(Bd4_A_East.marker)[1], function(x){
    if (Bd4_A_East.marker$AlleleCount1[x] == 2 | Bd4_A_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd4_A_East_merge <- merge(Bd4_A_East.marker, Bd4_A_East_Clock_J_age_estimate, by = "MarkerID")
  
  Bd5_A_East.marker$Frequencies <- Bd5_A_East.marker$AlleleCount1/(Bd5_A_East.marker$AlleleCount1 + Bd5_A_East.marker$AlleleCount0)
  Bd5_A_East.marker$Sigelton <- unlist(sapply(1:dim(Bd5_A_East.marker)[1], function(x){
    if (Bd5_A_East.marker$AlleleCount1[x] == 2 | Bd5_A_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd5_A_East_merge <- merge(Bd5_A_East.marker, Bd5_A_East_Clock_J_age_estimate, by = "MarkerID")
  
  All_A_East_merge <- rbind(Bd1_A_East_merge, Bd2_A_East_merge, Bd3_A_East_merge, Bd4_A_East_merge, Bd5_A_East_merge)
  
  synonymous_Bd1_A_East_position <- as.numeric(synonymous_A_East.vcf@fix[synonymous_A_East.vcf@fix[,1]=="Bd1",2])
  synonymous_Bd1_A_East_age <- Bd1_A_East_merge[Bd1_A_East_merge[,3]%in%synonymous_Bd1_A_East_position, c(1:3, 15, 20)]
  
  synonymous_Bd2_A_East_position <- as.numeric(synonymous_A_East.vcf@fix[synonymous_A_East.vcf@fix[,1]=="Bd2",2])
  synonymous_Bd2_A_East_age <- Bd2_A_East_merge[Bd2_A_East_merge[,3]%in%synonymous_Bd2_A_East_position, c(1:3, 15, 20)]
  
  synonymous_Bd3_A_East_position <- as.numeric(synonymous_A_East.vcf@fix[synonymous_A_East.vcf@fix[,1]=="Bd3",2])
  synonymous_Bd3_A_East_age <- Bd3_A_East_merge[Bd3_A_East_merge[,3]%in%synonymous_Bd3_A_East_position, c(1:3, 15, 20)]
  
  synonymous_Bd4_A_East_position <- as.numeric(synonymous_A_East.vcf@fix[synonymous_A_East.vcf@fix[,1]=="Bd4",2])
  synonymous_Bd4_A_East_age <- Bd4_A_East_merge[Bd4_A_East_merge[,3]%in%synonymous_Bd4_A_East_position, c(1:3, 15, 20)]
  
  synonymous_Bd5_A_East_position <- as.numeric(synonymous_A_East.vcf@fix[synonymous_A_East.vcf@fix[,1]=="Bd5",2])
  synonymous_Bd5_A_East_age <- Bd5_A_East_merge[Bd5_A_East_merge[,3]%in%synonymous_Bd5_A_East_position, c(1:3, 15, 20)]
  
  synonymous_All_A_East_age <- rbind(synonymous_Bd1_A_East_age, synonymous_Bd2_A_East_age, synonymous_Bd3_A_East_age, synonymous_Bd4_A_East_age, synonymous_Bd5_A_East_age)
  
}
## B West
{
  synonymous_B_West.vcf <- read.vcfR("./VCF_files/four_fold_degenerate_SNP_Bdis326_B_West.vcf", verbose = FALSE )
  
  Bd1_B_West_Clock_J_age_estimate <- read.table("./GEVA/B_West/Bd1_B_West_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd2_B_West_Clock_J_age_estimate <- read.table("./GEVA/B_West/Bd2_B_West_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd3_B_West_Clock_J_age_estimate <- read.table("./GEVA/B_West/Bd3_B_West_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd4_B_West_Clock_J_age_estimate <- read.table("./GEVA/B_West/Bd4_B_West_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd5_B_West_Clock_J_age_estimate <- read.table("./GEVA/B_West/Bd5_B_West_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd1_B_West.marker <- read.table("./GEVA/B_West/Bd1_B_West.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd2_B_West.marker <- read.table("./GEVA/B_West/Bd2_B_West.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd3_B_West.marker <- read.table("./GEVA/B_West/Bd3_B_West.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd4_B_West.marker <- read.table("./GEVA/B_West/Bd4_B_West.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd5_B_West.marker <- read.table("./GEVA/B_West/Bd5_B_West.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  not_derived_TE_Bd1_B_West <- read.table("./GEVA/B_West/not_derived_TE_Bd1_B_West.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd2_B_West <- read.table("./GEVA/B_West/not_derived_TE_Bd2_B_West.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd3_B_West <- read.table("./GEVA/B_West/not_derived_TE_Bd3_B_West.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd4_B_West <- read.table("./GEVA/B_West/not_derived_TE_Bd4_B_West.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd5_B_West <- read.table("./GEVA/B_West/not_derived_TE_Bd5_B_West.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  not_derived_TE_All_B_West <- rbind(not_derived_TE_Bd1_B_West, not_derived_TE_Bd2_B_West, not_derived_TE_Bd3_B_West, not_derived_TE_Bd4_B_West, not_derived_TE_Bd5_B_West)  
  
  Bd1_B_West.marker$Frequencies <- Bd1_B_West.marker$AlleleCount1/(Bd1_B_West.marker$AlleleCount1 + Bd1_B_West.marker$AlleleCount0)
  Bd1_B_West.marker$Sigelton <- unlist(sapply(1:dim(Bd1_B_West.marker)[1], function(x){
    if (Bd1_B_West.marker$AlleleCount1[x] == 2 | Bd1_B_West.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd1_B_West_merge <- merge(Bd1_B_West.marker, Bd1_B_West_Clock_J_age_estimate, by = "MarkerID")
  
  Bd2_B_West.marker$Frequencies <- Bd2_B_West.marker$AlleleCount1/(Bd2_B_West.marker$AlleleCount1 + Bd2_B_West.marker$AlleleCount0)
  Bd2_B_West.marker$Sigelton <- unlist(sapply(1:dim(Bd2_B_West.marker)[1], function(x){
    if (Bd2_B_West.marker$AlleleCount1[x] == 2 | Bd2_B_West.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd2_B_West_merge <- merge(Bd2_B_West.marker, Bd2_B_West_Clock_J_age_estimate, by = "MarkerID")
  
  Bd3_B_West.marker$Frequencies <- Bd3_B_West.marker$AlleleCount1/(Bd3_B_West.marker$AlleleCount1 + Bd3_B_West.marker$AlleleCount0)
  Bd3_B_West.marker$Sigelton <- unlist(sapply(1:dim(Bd3_B_West.marker)[1], function(x){
    if (Bd3_B_West.marker$AlleleCount1[x] == 2 | Bd3_B_West.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd3_B_West_merge <- merge(Bd3_B_West.marker, Bd3_B_West_Clock_J_age_estimate, by = "MarkerID")
  
  Bd4_B_West.marker$Frequencies <- Bd4_B_West.marker$AlleleCount1/(Bd4_B_West.marker$AlleleCount1 + Bd4_B_West.marker$AlleleCount0)
  Bd4_B_West.marker$Sigelton <- unlist(sapply(1:dim(Bd4_B_West.marker)[1], function(x){
    if (Bd4_B_West.marker$AlleleCount1[x] == 2 | Bd4_B_West.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd4_B_West_merge <- merge(Bd4_B_West.marker, Bd4_B_West_Clock_J_age_estimate, by = "MarkerID")
  
  Bd5_B_West.marker$Frequencies <- Bd5_B_West.marker$AlleleCount1/(Bd5_B_West.marker$AlleleCount1 + Bd5_B_West.marker$AlleleCount0)
  Bd5_B_West.marker$Sigelton <- unlist(sapply(1:dim(Bd5_B_West.marker)[1], function(x){
    if (Bd5_B_West.marker$AlleleCount1[x] == 2 | Bd5_B_West.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd5_B_West_merge <- merge(Bd5_B_West.marker, Bd5_B_West_Clock_J_age_estimate, by = "MarkerID")
  
  All_B_West_merge <- rbind(Bd1_B_West_merge, Bd2_B_West_merge, Bd3_B_West_merge, Bd4_B_West_merge, Bd5_B_West_merge)
  
  synonymous_Bd1_B_West_position <- as.numeric(synonymous_B_West.vcf@fix[synonymous_B_West.vcf@fix[,1]=="Bd1",2])
  synonymous_Bd1_B_West_age <- Bd1_B_West_merge[Bd1_B_West_merge[,3]%in%synonymous_Bd1_B_West_position, c(1:3, 15, 20)]
  
  synonymous_Bd2_B_West_position <- as.numeric(synonymous_B_West.vcf@fix[synonymous_B_West.vcf@fix[,1]=="Bd2",2])
  synonymous_Bd2_B_West_age <- Bd2_B_West_merge[Bd2_B_West_merge[,3]%in%synonymous_Bd2_B_West_position, c(1:3, 15, 20)]
  
  synonymous_Bd3_B_West_position <- as.numeric(synonymous_B_West.vcf@fix[synonymous_B_West.vcf@fix[,1]=="Bd3",2])
  synonymous_Bd3_B_West_age <- Bd3_B_West_merge[Bd3_B_West_merge[,3]%in%synonymous_Bd3_B_West_position, c(1:3, 15, 20)]
  
  synonymous_Bd4_B_West_position <- as.numeric(synonymous_B_West.vcf@fix[synonymous_B_West.vcf@fix[,1]=="Bd4",2])
  synonymous_Bd4_B_West_age <- Bd4_B_West_merge[Bd4_B_West_merge[,3]%in%synonymous_Bd4_B_West_position, c(1:3, 15, 20)]
  
  synonymous_Bd5_B_West_position <- as.numeric(synonymous_B_West.vcf@fix[synonymous_B_West.vcf@fix[,1]=="Bd5",2])
  synonymous_Bd5_B_West_age <- Bd5_B_West_merge[Bd5_B_West_merge[,3]%in%synonymous_Bd5_B_West_position, c(1:3, 15, 20)]
  
  synonymous_All_B_West_age <- rbind(synonymous_Bd1_B_West_age, synonymous_Bd2_B_West_age, synonymous_Bd3_B_West_age, synonymous_Bd4_B_West_age, synonymous_Bd5_B_West_age)
  
}
## B East
{
  synonymous_B_East.vcf <- read.vcfR("./VCF_files/four_fold_degenerate_SNP_Bdis326_B_East.vcf", verbose = FALSE )
  
  Bd1_B_East_Clock_J_age_estimate <- read.table("./GEVA/B_East/Bd1_B_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd2_B_East_Clock_J_age_estimate <- read.table("./GEVA/B_East/Bd2_B_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd3_B_East_Clock_J_age_estimate <- read.table("./GEVA/B_East/Bd3_B_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd4_B_East_Clock_J_age_estimate <- read.table("./GEVA/B_East/Bd4_B_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd5_B_East_Clock_J_age_estimate <- read.table("./GEVA/B_East/Bd5_B_East_Clock_J_age_estimate.sites2.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd1_B_East.marker <- read.table("./GEVA/B_East/Bd1_B_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd2_B_East.marker <- read.table("./GEVA/B_East/Bd2_B_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd3_B_East.marker <- read.table("./GEVA/B_East/Bd3_B_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd4_B_East.marker <- read.table("./GEVA/B_East/Bd4_B_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  Bd5_B_East.marker <- read.table("./GEVA/B_East/Bd5_B_East.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  not_derived_TE_Bd1_B_East <- read.table("./GEVA/B_East/not_derived_TE_Bd1_B_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd2_B_East <- read.table("./GEVA/B_East/not_derived_TE_Bd2_B_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd3_B_East <- read.table("./GEVA/B_East/not_derived_TE_Bd3_B_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd4_B_East <- read.table("./GEVA/B_East/not_derived_TE_Bd4_B_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  not_derived_TE_Bd5_B_East <- read.table("./GEVA/B_East/not_derived_TE_Bd5_B_East.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  not_derived_TE_All_B_East <- rbind(not_derived_TE_Bd1_B_East, not_derived_TE_Bd2_B_East, not_derived_TE_Bd3_B_East, not_derived_TE_Bd4_B_East, not_derived_TE_Bd5_B_East)  
  
  Bd1_B_East.marker$Frequencies <- Bd1_B_East.marker$AlleleCount1/(Bd1_B_East.marker$AlleleCount1 + Bd1_B_East.marker$AlleleCount0)
  Bd1_B_East.marker$Sigelton <- unlist(sapply(1:dim(Bd1_B_East.marker)[1], function(x){
    if (Bd1_B_East.marker$AlleleCount1[x] == 2 | Bd1_B_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd1_B_East_merge <- merge(Bd1_B_East.marker, Bd1_B_East_Clock_J_age_estimate, by = "MarkerID")
  
  Bd2_B_East.marker$Frequencies <- Bd2_B_East.marker$AlleleCount1/(Bd2_B_East.marker$AlleleCount1 + Bd2_B_East.marker$AlleleCount0)
  Bd2_B_East.marker$Sigelton <- unlist(sapply(1:dim(Bd2_B_East.marker)[1], function(x){
    if (Bd2_B_East.marker$AlleleCount1[x] == 2 | Bd2_B_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd2_B_East_merge <- merge(Bd2_B_East.marker, Bd2_B_East_Clock_J_age_estimate, by = "MarkerID")
  
  Bd3_B_East.marker$Frequencies <- Bd3_B_East.marker$AlleleCount1/(Bd3_B_East.marker$AlleleCount1 + Bd3_B_East.marker$AlleleCount0)
  Bd3_B_East.marker$Sigelton <- unlist(sapply(1:dim(Bd3_B_East.marker)[1], function(x){
    if (Bd3_B_East.marker$AlleleCount1[x] == 2 | Bd3_B_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd3_B_East_merge <- merge(Bd3_B_East.marker, Bd3_B_East_Clock_J_age_estimate, by = "MarkerID")
  
  Bd4_B_East.marker$Frequencies <- Bd4_B_East.marker$AlleleCount1/(Bd4_B_East.marker$AlleleCount1 + Bd4_B_East.marker$AlleleCount0)
  Bd4_B_East.marker$Sigelton <- unlist(sapply(1:dim(Bd4_B_East.marker)[1], function(x){
    if (Bd4_B_East.marker$AlleleCount1[x] == 2 | Bd4_B_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd4_B_East_merge <- merge(Bd4_B_East.marker, Bd4_B_East_Clock_J_age_estimate, by = "MarkerID")
  
  Bd5_B_East.marker$Frequencies <- Bd5_B_East.marker$AlleleCount1/(Bd5_B_East.marker$AlleleCount1 + Bd5_B_East.marker$AlleleCount0)
  Bd5_B_East.marker$Sigelton <- unlist(sapply(1:dim(Bd5_B_East.marker)[1], function(x){
    if (Bd5_B_East.marker$AlleleCount1[x] == 2 | Bd5_B_East.marker$AlleleCount0[x] == 2 ) {
      return("Yes")
    } else {return("No")}
  }))
  Bd5_B_East_merge <- merge(Bd5_B_East.marker, Bd5_B_East_Clock_J_age_estimate, by = "MarkerID")
  
  All_B_East_merge <- rbind(Bd1_B_East_merge, Bd2_B_East_merge, Bd3_B_East_merge, Bd4_B_East_merge, Bd5_B_East_merge)
  
  synonymous_Bd1_B_East_position <- as.numeric(synonymous_B_East.vcf@fix[synonymous_B_East.vcf@fix[,1]=="Bd1",2])
  synonymous_Bd1_B_East_age <- Bd1_B_East_merge[Bd1_B_East_merge[,3]%in%synonymous_Bd1_B_East_position, c(1:3, 15, 20)]
  
  synonymous_Bd2_B_East_position <- as.numeric(synonymous_B_East.vcf@fix[synonymous_B_East.vcf@fix[,1]=="Bd2",2])
  synonymous_Bd2_B_East_age <- Bd2_B_East_merge[Bd2_B_East_merge[,3]%in%synonymous_Bd2_B_East_position, c(1:3, 15, 20)]
  
  synonymous_Bd3_B_East_position <- as.numeric(synonymous_B_East.vcf@fix[synonymous_B_East.vcf@fix[,1]=="Bd3",2])
  synonymous_Bd3_B_East_age <- Bd3_B_East_merge[Bd3_B_East_merge[,3]%in%synonymous_Bd3_B_East_position, c(1:3, 15, 20)]
  
  synonymous_Bd4_B_East_position <- as.numeric(synonymous_B_East.vcf@fix[synonymous_B_East.vcf@fix[,1]=="Bd4",2])
  synonymous_Bd4_B_East_age <- Bd4_B_East_merge[Bd4_B_East_merge[,3]%in%synonymous_Bd4_B_East_position, c(1:3, 15, 20)]
  
  synonymous_Bd5_B_East_position <- as.numeric(synonymous_B_East.vcf@fix[synonymous_B_East.vcf@fix[,1]=="Bd5",2])
  synonymous_Bd5_B_East_age <- Bd5_B_East_merge[Bd5_B_East_merge[,3]%in%synonymous_Bd5_B_East_position, c(1:3, 15, 20)]
  
  synonymous_All_B_East_age <- rbind(synonymous_Bd1_B_East_age, synonymous_Bd2_B_East_age, synonymous_Bd3_B_East_age, synonymous_Bd4_B_East_age, synonymous_Bd5_B_East_age)
}

## define some functions
downsampling.neutral.SNP.function <- function(SNP.data, TE.data){
  resampled_SNP_age <- matrix(NA, nrow = dim(TE.data)[1], ncol = dim(SNP.data)[2])
  colnames(resampled_SNP_age) <- colnames(SNP.data)
  
  
  for (s in 1:dim(TE.data)[1]) {
    my_age <- TE.data[s,]$PostMode
    same_age_tab_b <- SNP.data[SNP.data$PostMode <= my_age + 100 & SNP.data$PostMode >= my_age - 100,]
    same_age_tab <- same_age_tab_b[!same_age_tab_b[,1]%in%resampled_SNP_age[,1],]
    
    
    if (dim(same_age_tab)[1] > 1) {
      resampled_SNP_age[s,] <- unlist(same_age_tab[sample(1:dim(same_age_tab)[1], 1),])
    } else if (dim(same_age_tab)[1] == 1) {
      resampled_SNP_age[s,] <- unlist(same_age_tab)
    } else {
      my_loop <- 2
      while(my_loop < 1000 & is.na(resampled_SNP_age[s,1])){
        
        same_age_tab_b <- SNP.data[SNP.data$PostMode <= my_age + 100*my_loop & SNP.data$PostMode >= my_age - 100*my_loop,]
        same_age_tab <- same_age_tab_b[!same_age_tab_b[,1]%in%resampled_SNP_age[,1],]
        
        if (dim(same_age_tab)[1] > 1) {
          resampled_SNP_age[s,] <- unlist(same_age_tab[sample(1:dim(same_age_tab)[1], 1),])
        } else if (dim(same_age_tab)[1] == 1) {
          resampled_SNP_age[s,] <- unlist(same_age_tab)
        }
        
        my_loop <- my_loop + 1
      }
    }
  }
  
  return(resampled_SNP_age)
  
  
}
decile.delta.age.fix.decile.size.function <- function(resampled.SNP, TE.data, decile.age.bonderies){
  
  Bonderies <- c(decile.age.bonderies[1:10], decile.age.bonderies[11]+1 )
  
  delta_freq <- c()
  for (s in 1:10) {
    delta_freq <- c(delta_freq, 
                    mean(TE.data[TE.data$PostMode >= Bonderies[s] & TE.data$PostMode < Bonderies[s+1], ]$Frequencies) - 
                      mean(resampled.SNP[resampled.SNP[,5] >= Bonderies[s] & resampled.SNP[,5] < Bonderies[s+1], 4])
    )
  }
  
  return(delta_freq)
}
decile.delta.age.fix.decile.size.variable.n.bin.function <- function(resampled.SNP, TE.data, decile.age.bonderies){
  n_bin <- length(decile.age.bonderies)
  Bonderies <- c(decile.age.bonderies[1:(n_bin-1)], decile.age.bonderies[n_bin]+1 )
  
  delta_freq <- c()
  for (s in 1:(n_bin-1)) {
    delta_freq <- c(delta_freq, 
                    mean(TE.data[TE.data$PostMode >= Bonderies[s] & TE.data$PostMode < Bonderies[s+1], ]$Frequencies) - 
                      mean(resampled.SNP[resampled.SNP[,5] >= Bonderies[s] & resampled.SNP[,5] < Bonderies[s+1], 4], na.rm = TRUE)
    )
  }
  
  return(delta_freq)
}
analyse.delta.age.bootstrapp.three.TEs.groups.all.genome.function <- function(n_boot, Clade, TE_ID_list, TE_ID_list_2, TE_ID_list_3, Name_list, Name_list_2, Name_list_3, Plot_Name, output_path = "/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/Age_analyses_results/", Nbin = 10, Singelton = FALSE){
  my_base_names <- c("All_CLADE_merge",
                     "not_derived_TE_All_CLADE",
                     "synonymous_All_CLADE_age")
  my_name_list <- gsub(pattern = "CLADE", replacement = Clade, x = my_base_names)
  
  ## get estimates
  if(Singelton){
    Retrotransposons_All_age_A <- get(my_name_list[1])[get(my_name_list[1])$Label %in% TE_ID_list & !get(my_name_list[1])$Position %in% get(my_name_list[2])[,2] & get(my_name_list[1])$Sigelton == "No",c(1:4,15,20)]
    Retrotransposons_All_age_B <- get(my_name_list[1])[get(my_name_list[1])$Label %in% TE_ID_list_2 & !get(my_name_list[1])$Position %in% get(my_name_list[2])[,2] & get(my_name_list[1])$Sigelton == "No",c(1:4,15,20)]
    Retrotransposons_All_age_C <- get(my_name_list[1])[get(my_name_list[1])$Label %in% TE_ID_list_3 & !get(my_name_list[1])$Position %in% get(my_name_list[2])[,2] & get(my_name_list[1])$Sigelton == "No",c(1:4,15,20)]
  } else {
    Retrotransposons_All_age_A <- get(my_name_list[1])[get(my_name_list[1])$Label %in% TE_ID_list & !get(my_name_list[1])$Position %in% get(my_name_list[2])[,2],c(1:4,15,20)]
    Retrotransposons_All_age_B <- get(my_name_list[1])[get(my_name_list[1])$Label %in% TE_ID_list_2 & !get(my_name_list[1])$Position %in% get(my_name_list[2])[,2],c(1:4,15,20)]
    Retrotransposons_All_age_C <- get(my_name_list[1])[get(my_name_list[1])$Label %in% TE_ID_list_3 & !get(my_name_list[1])$Position %in% get(my_name_list[2])[,2],c(1:4,15,20)]
  }
  
  # downsapled neutral data 
  resampled_synonymous_All_age_A <- downsampling.neutral.SNP.function(get(my_name_list[3]), Retrotransposons_All_age_A)
  resampled_synonymous_All_age_B <- downsampling.neutral.SNP.function(get(my_name_list[3]), Retrotransposons_All_age_B)
  resampled_synonymous_All_age_C <- downsampling.neutral.SNP.function(get(my_name_list[3]), Retrotransposons_All_age_C)
  
  
  ## sort data 
  Retrotransposons_All_age_sorted_A <- Retrotransposons_All_age_A[order(Retrotransposons_All_age_A$PostMode), ]
  resampled_synonymous_All_age_sorted_A <- resampled_synonymous_All_age_A[order(resampled_synonymous_All_age_A[,5]),]
  
  Retrotransposons_All_age_sorted_B <- Retrotransposons_All_age_B[order(Retrotransposons_All_age_B$PostMode), ]
  resampled_synonymous_All_age_sorted_B <- resampled_synonymous_All_age_B[order(resampled_synonymous_All_age_B[,5]),]
  
  Retrotransposons_All_age_sorted_C <- Retrotransposons_All_age_C[order(Retrotransposons_All_age_C$PostMode), ]
  resampled_synonymous_All_age_sorted_C <- resampled_synonymous_All_age_C[order(resampled_synonymous_All_age_C[,5]),]
  
  Retrotransposons_All_age_Full <- rbind(Retrotransposons_All_age_A, Retrotransposons_All_age_B, Retrotransposons_All_age_C)
  resampled_synonymous_All_age_Full <- rbind(resampled_synonymous_All_age_A, resampled_synonymous_All_age_B, resampled_synonymous_All_age_C)
  Retrotransposons_All_age_sorted_Full <- Retrotransposons_All_age_Full[order(Retrotransposons_All_age_Full$PostMode), ]
  resampled_synonymous_All_age_sorted_Full <- resampled_synonymous_All_age_Full[order(resampled_synonymous_All_age_Full[,5]),]
  
  
  ## get delta age:
  All_ages <- sort(c(Retrotransposons_All_age_sorted_Full$PostMode, resampled_synonymous_All_age_sorted_Full[,5]))
  my_age_n <- length(All_ages)
  my_bonderies <- c(0)
  for (n in 1:Nbin) {
    my_bonderies <- c(my_bonderies, All_ages[floor(n*my_age_n/Nbin)])
    
  }
  
  ## bootstrapping
  delta_freq_decile_All_boot_A <- matrix(unlist(sapply(1:n_boot, function(x){
    Retrotransposons_All_age_resample <- Retrotransposons_All_age_A[sample(1:dim(Retrotransposons_All_age_A)[1], dim(Retrotransposons_All_age_A)[1], replace = TRUE),]
    resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(get(my_name_list[3]), Retrotransposons_All_age_resample)
    Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
    resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
    delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies)
    return(delta_freq_decile_All_resample)
  }))
  , nrow = n_boot, ncol = Nbin, byrow = TRUE)
  
  delta_freq_decile_All_boot_B <- matrix(unlist(sapply(1:n_boot, function(x){
    Retrotransposons_All_age_resample <- Retrotransposons_All_age_B[sample(1:dim(Retrotransposons_All_age_B)[1], dim(Retrotransposons_All_age_B)[1], replace = TRUE),]
    resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(get(my_name_list[3]), Retrotransposons_All_age_resample)
    Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
    resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
    delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies)
    return(delta_freq_decile_All_resample)
  }))
  , nrow = n_boot, ncol = Nbin, byrow = TRUE)
  
  delta_freq_decile_All_boot_C <- matrix(unlist(sapply(1:n_boot, function(x){
    Retrotransposons_All_age_resample <- Retrotransposons_All_age_C[sample(1:dim(Retrotransposons_All_age_C)[1], dim(Retrotransposons_All_age_C)[1], replace = TRUE),]
    resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(get(my_name_list[3]), Retrotransposons_All_age_resample)
    Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
    resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
    delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies)
    return(delta_freq_decile_All_resample)
  }))
  , nrow = n_boot, ncol = Nbin, byrow = TRUE)
  
  my_plot_seq <- as.vector(rbind(c(1:Nbin), c((Nbin+1):(2*Nbin)), c((2*Nbin+1):(3*Nbin))))
  my_at <- as.vector(rbind(seq(1.2, 3*Nbin, by = 3), 
                           seq(2, 3*Nbin, by = 3), 
                           seq(2.8, 3*Nbin, by = 3)))
  my_labels <- c()
  for (h in 1:Nbin) {
    my_labels <- c(my_labels, paste0(round(my_bonderies[h], 0), "\n-\n", round(my_bonderies[h+1], 0)))
  }
  
  ## plotting
  pdf(file = paste0(output_path, "Plot_delta_age_deciles_", Clade, "_", Plot_Name, "_", Nbin, "_", n_boot, "_bootstrapping.pdf"), width = 11)
  
  boxplot(cbind(delta_freq_decile_All_boot_A, delta_freq_decile_All_boot_B, delta_freq_decile_All_boot_C)[,my_plot_seq],
          col = c(2,3,4), xlab = "Age bins", ylab = "delta frequency", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, main = paste0("Retrotransposons ", gsub(pattern = "_", replacement = " ", Clade)), 
          names = c(1:Nbin,1:Nbin,1:Nbin), xaxt = "n", 
          at = my_at)
  abline(h = 0, lty = 3)
  axis(1, at=seq(2,3*Nbin,3), line = 1, tick = FALSE, labels = my_labels)
  
  
  legend("bottomleft", legend = c(Name_list, Name_list_2, Name_list_3), col = c(2,3,4), pch = 15, bty = "n", cex = 2)
  dev.off()
  
}
match.SNP.frequency.to.TE <- function(SNP.data, TE.data, bin.size = 0.01){
  my_breaks <- seq(0,1,bin.size)
  my_TE_counts <- unlist(sapply(1:100, function(x){
    return(dim(TE.data[TE.data$Frequencies < my_breaks[x+1] & TE.data$Frequencies >= my_breaks[x], ])[1])
  }))
  my_SNP_counts <- unlist(sapply(1:100, function(x){
    return(dim(SNP.data[SNP.data$Frequencies < my_breaks[x+1] & SNP.data$Frequencies >= my_breaks[x], ])[1])
  }))
  my_prop <- my_TE_counts/sum(my_TE_counts)
  my_exp_SNP_numbers <- ceiling(sum(my_SNP_counts)*my_prop)
  my_dif <- (my_SNP_counts - my_exp_SNP_numbers)/my_SNP_counts
  my_lim <- max(1:100*(my_dif==min(my_dif, na.rm = TRUE)), na.rm = TRUE)
  my_target_SNP_numbers <- ceiling(floor(my_SNP_counts[my_lim]/my_prop[my_lim])*my_prop)
  
  
  my_frequency_matched_SNP_ID <- unlist(sapply(1:100, function(x){
    my_full_data <- SNP.data[SNP.data$Frequencies < my_breaks[x+1] & SNP.data$Frequencies >= my_breaks[x], ]
    my_freq_mat <- my_full_data[sample(1:dim(my_full_data)[1], size = my_target_SNP_numbers[x], replace = FALSE),]
    return(my_freq_mat$MarkerID)
  }))
  
  my_frequency_matched_SNP <- SNP.data[SNP.data$MarkerID %in% my_frequency_matched_SNP_ID,]
  
  return(my_frequency_matched_SNP)
}


############################################
######### the analyses   ###################
############################################

## get the different list of TEs:
Retrotransposons_All_age_B_East <- All_B_East_merge[All_B_East_merge$Label %in% Retrotransposon_ID & !All_B_East_merge$Position %in% not_derived_TE_All_B_East[,2],c(1:4,15,20)]
Retrotransposons_All_age_A_Italia <- All_A_Italia_merge[All_A_Italia_merge$Label %in% Retrotransposon_ID & !All_A_Italia_merge$Position %in% not_derived_TE_All_A_Italia[,2],c(1:4,15,20)]
Retrotransposons_All_age_A_East <- All_A_East_merge[All_A_East_merge$Label %in% Retrotransposon_ID & !All_A_East_merge$Position %in% not_derived_TE_All_A_East[,2],c(1:4,15,20)]
Retrotransposons_All_age_B_West <- All_B_West_merge[All_B_West_merge$Label %in% Retrotransposon_ID & !All_B_West_merge$Position %in% not_derived_TE_All_B_West[,2],c(1:4,15,20)]

## match frequency distribution:
freq_matched_synonymous_All_age_B_East <- match.SNP.frequency.to.TE(synonymous_All_B_East_age, Retrotransposons_All_age_B_East)
freq_matched_synonymous_All_age_A_Italia <- match.SNP.frequency.to.TE(synonymous_All_A_Italia_age, Retrotransposons_All_age_A_Italia)
freq_matched_synonymous_All_age_A_East <- match.SNP.frequency.to.TE(synonymous_All_A_East_age, Retrotransposons_All_age_A_East)
freq_matched_synonymous_All_age_B_West <- match.SNP.frequency.to.TE(synonymous_All_B_West_age, Retrotransposons_All_age_B_West)

# downsapled neutral data 
resampled_freq_matched_synonymous_All_age_B_East <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_East, Retrotransposons_All_age_B_East)
resampled_freq_matched_synonymous_All_age_A_Italia <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_Italia, Retrotransposons_All_age_A_Italia)
resampled_freq_matched_synonymous_All_age_A_East <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_East, Retrotransposons_All_age_A_East)
resampled_freq_matched_synonymous_All_age_B_West <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_West, Retrotransposons_All_age_B_West)

## sort data 
Retrotransposons_All_age_sorted_B_East <- Retrotransposons_All_age_B_East[order(Retrotransposons_All_age_B_East$PostMode), ]
resampled_freq_matched_synonymous_All_age_sorted_B_East <- resampled_freq_matched_synonymous_All_age_B_East[order(resampled_freq_matched_synonymous_All_age_B_East[,5]),]
Retrotransposons_All_age_sorted_A_Italia <- Retrotransposons_All_age_A_Italia[order(Retrotransposons_All_age_A_Italia$PostMode), ]
resampled_freq_matched_synonymous_All_age_sorted_A_Italia <- resampled_freq_matched_synonymous_All_age_A_Italia[order(resampled_freq_matched_synonymous_All_age_A_Italia[,5]),]
Retrotransposons_All_age_sorted_A_East <- Retrotransposons_All_age_A_East[order(Retrotransposons_All_age_A_East$PostMode), ]
resampled_freq_matched_synonymous_All_age_sorted_A_East <- resampled_freq_matched_synonymous_All_age_A_East[order(resampled_freq_matched_synonymous_All_age_A_East[,5]),]
Retrotransposons_All_age_sorted_B_West <- Retrotransposons_All_age_B_West[order(Retrotransposons_All_age_B_West$PostMode), ]
resampled_freq_matched_synonymous_All_age_sorted_B_West <- resampled_freq_matched_synonymous_All_age_B_West[order(resampled_freq_matched_synonymous_All_age_B_West[,5]),]

## get bins
Nbin = 15
All_ages_B_East <- sort(c(Retrotransposons_All_age_sorted_B_East$PostMode, resampled_freq_matched_synonymous_All_age_sorted_B_East[,5]))
my_age_n_B_East <- length(All_ages_B_East)
my_bonderies_B_East <- c(0)
for (n in 1:Nbin) {
  my_bonderies_B_East <- c(my_bonderies_B_East, All_ages_B_East[floor(n*my_age_n_B_East/Nbin)])
  
}
All_ages_A_Italia <- sort(c(Retrotransposons_All_age_sorted_A_Italia$PostMode, resampled_freq_matched_synonymous_All_age_sorted_A_Italia[,5]))
my_age_n_A_Italia <- length(All_ages_A_Italia)
my_bonderies_A_Italia <- c(0)
for (n in 1:Nbin) {
  my_bonderies_A_Italia <- c(my_bonderies_A_Italia, All_ages_A_Italia[floor(n*my_age_n_A_Italia/Nbin)])
  
}
All_ages_A_East <- sort(c(Retrotransposons_All_age_sorted_A_East$PostMode, resampled_freq_matched_synonymous_All_age_sorted_A_East[,5]))
my_age_n_A_East <- length(All_ages_A_East)
my_bonderies_A_East <- c(0)
for (n in 1:Nbin) {
  my_bonderies_A_East <- c(my_bonderies_A_East, All_ages_A_East[floor(n*my_age_n_A_East/Nbin)])
  
}
All_ages_B_West <- sort(c(Retrotransposons_All_age_sorted_B_West$PostMode, resampled_freq_matched_synonymous_All_age_sorted_B_West[,5]))
my_age_n_B_West <- length(All_ages_B_West)
my_bonderies_B_West <- c(0)
for (n in 1:Nbin) {
  my_bonderies_B_West <- c(my_bonderies_B_West, All_ages_B_West[floor(n*my_age_n_B_West/Nbin)])
  
}

## get delta frequency:
delta_freq_decile_All_resample_freq_matched_B_East <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_freq_matched_synonymous_All_age_sorted_B_East, Retrotransposons_All_age_sorted_B_East, my_bonderies_B_East)
delta_freq_decile_All_resample_freq_matched_A_Italia <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_freq_matched_synonymous_All_age_sorted_A_Italia, Retrotransposons_All_age_sorted_A_Italia, my_bonderies_A_Italia)
delta_freq_decile_All_resample_freq_matched_A_East <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_freq_matched_synonymous_All_age_sorted_A_East, Retrotransposons_All_age_sorted_A_East, my_bonderies_A_East)
delta_freq_decile_All_resample_freq_matched_B_West <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_freq_matched_synonymous_All_age_sorted_B_West, Retrotransposons_All_age_sorted_B_West, my_bonderies_B_West)

## bootstrapping
n_boot = 100
delta_freq_decile_All_boot_freq_matched_B_East <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_B_East[sample(1:dim(Retrotransposons_All_age_B_East)[1], dim(Retrotransposons_All_age_B_East)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_East, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_East)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_A_Italia <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_A_Italia[sample(1:dim(Retrotransposons_All_age_A_Italia)[1], dim(Retrotransposons_All_age_A_Italia)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_Italia, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_Italia)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_A_East <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_A_East[sample(1:dim(Retrotransposons_All_age_A_East)[1], dim(Retrotransposons_All_age_A_East)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_East, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_East)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_B_West <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_B_West[sample(1:dim(Retrotransposons_All_age_B_West)[1], dim(Retrotransposons_All_age_B_West)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_West, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_West)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)


## Analyses the snpeff results:

## read in data:
## B West
{
  high_effect_SNPs_B_West.vcf <- read.vcfR("./Bdis326_SNPs_filtered_snpEff_high_effect.vcf", verbose = FALSE )
  
  high_effect_SNPs_Bd1_B_West_position <- as.numeric(high_effect_SNPs_B_West.vcf@fix[high_effect_SNPs_B_West.vcf@fix[,1]=="1",2])
  high_effect_SNPs_Bd1_B_West_age <- Bd1_B_West_merge[Bd1_B_West_merge[,3]%in%high_effect_SNPs_Bd1_B_West_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd2_B_West_position <- as.numeric(high_effect_SNPs_B_West.vcf@fix[high_effect_SNPs_B_West.vcf@fix[,1]=="2",2])
  high_effect_SNPs_Bd2_B_West_age <- Bd2_B_West_merge[Bd2_B_West_merge[,3]%in%high_effect_SNPs_Bd2_B_West_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd3_B_West_position <- as.numeric(high_effect_SNPs_B_West.vcf@fix[high_effect_SNPs_B_West.vcf@fix[,1]=="3",2])
  high_effect_SNPs_Bd3_B_West_age <- Bd3_B_West_merge[Bd3_B_West_merge[,3]%in%high_effect_SNPs_Bd3_B_West_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd4_B_West_position <- as.numeric(high_effect_SNPs_B_West.vcf@fix[high_effect_SNPs_B_West.vcf@fix[,1]=="4",2])
  high_effect_SNPs_Bd4_B_West_age <- Bd4_B_West_merge[Bd4_B_West_merge[,3]%in%high_effect_SNPs_Bd4_B_West_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd5_B_West_position <- as.numeric(high_effect_SNPs_B_West.vcf@fix[high_effect_SNPs_B_West.vcf@fix[,1]=="5",2])
  high_effect_SNPs_Bd5_B_West_age <- Bd5_B_West_merge[Bd5_B_West_merge[,3]%in%high_effect_SNPs_Bd5_B_West_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_All_B_West_age_full <- as.data.frame(rbind(high_effect_SNPs_Bd1_B_West_age, high_effect_SNPs_Bd2_B_West_age, high_effect_SNPs_Bd3_B_West_age, high_effect_SNPs_Bd4_B_West_age, high_effect_SNPs_Bd5_B_West_age))
  
  
  non_synonymous_B_West.vcf <- read.vcfR("./zero_fold_degenerate_SNP_Bdis326_B_West.vcf", verbose = FALSE )
  
  non_synonymous_Bd1_B_West_position <- as.numeric(non_synonymous_B_West.vcf@fix[non_synonymous_B_West.vcf@fix[,1]=="Bd1",2])
  non_synonymous_Bd1_B_West_age <- Bd1_B_West_merge[Bd1_B_West_merge[,3]%in%non_synonymous_Bd1_B_West_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd2_B_West_position <- as.numeric(non_synonymous_B_West.vcf@fix[non_synonymous_B_West.vcf@fix[,1]=="Bd2",2])
  non_synonymous_Bd2_B_West_age <- Bd2_B_West_merge[Bd2_B_West_merge[,3]%in%non_synonymous_Bd2_B_West_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd3_B_West_position <- as.numeric(non_synonymous_B_West.vcf@fix[non_synonymous_B_West.vcf@fix[,1]=="Bd3",2])
  non_synonymous_Bd3_B_West_age <- Bd3_B_West_merge[Bd3_B_West_merge[,3]%in%non_synonymous_Bd3_B_West_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd4_B_West_position <- as.numeric(non_synonymous_B_West.vcf@fix[non_synonymous_B_West.vcf@fix[,1]=="Bd4",2])
  non_synonymous_Bd4_B_West_age <- Bd4_B_West_merge[Bd4_B_West_merge[,3]%in%non_synonymous_Bd4_B_West_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd5_B_West_position <- as.numeric(non_synonymous_B_West.vcf@fix[non_synonymous_B_West.vcf@fix[,1]=="Bd5",2])
  non_synonymous_Bd5_B_West_age <- Bd5_B_West_merge[Bd5_B_West_merge[,3]%in%non_synonymous_Bd5_B_West_position, c(1:3, 15, 20)]
  
  non_synonymous_All_B_West_age_full <- as.data.frame(rbind(non_synonymous_Bd1_B_West_age, non_synonymous_Bd2_B_West_age, non_synonymous_Bd3_B_West_age, non_synonymous_Bd4_B_West_age, non_synonymous_Bd5_B_West_age))
  
}
## A Italia
{
  high_effect_SNPs_A_Italia.vcf <- read.vcfR("./Bdis326_SNPs_filtered_snpEff_high_effect.vcf", verbose = FALSE )
  
  high_effect_SNPs_Bd1_A_Italia_position <- as.numeric(high_effect_SNPs_A_Italia.vcf@fix[high_effect_SNPs_A_Italia.vcf@fix[,1]=="1",2])
  high_effect_SNPs_Bd1_A_Italia_age <- Bd1_A_Italia_merge[Bd1_A_Italia_merge[,3]%in%high_effect_SNPs_Bd1_A_Italia_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd2_A_Italia_position <- as.numeric(high_effect_SNPs_A_Italia.vcf@fix[high_effect_SNPs_A_Italia.vcf@fix[,1]=="2",2])
  high_effect_SNPs_Bd2_A_Italia_age <- Bd2_A_Italia_merge[Bd2_A_Italia_merge[,3]%in%high_effect_SNPs_Bd2_A_Italia_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd3_A_Italia_position <- as.numeric(high_effect_SNPs_A_Italia.vcf@fix[high_effect_SNPs_A_Italia.vcf@fix[,1]=="3",2])
  high_effect_SNPs_Bd3_A_Italia_age <- Bd3_A_Italia_merge[Bd3_A_Italia_merge[,3]%in%high_effect_SNPs_Bd3_A_Italia_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd4_A_Italia_position <- as.numeric(high_effect_SNPs_A_Italia.vcf@fix[high_effect_SNPs_A_Italia.vcf@fix[,1]=="4",2])
  high_effect_SNPs_Bd4_A_Italia_age <- Bd4_A_Italia_merge[Bd4_A_Italia_merge[,3]%in%high_effect_SNPs_Bd4_A_Italia_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd5_A_Italia_position <- as.numeric(high_effect_SNPs_A_Italia.vcf@fix[high_effect_SNPs_A_Italia.vcf@fix[,1]=="5",2])
  high_effect_SNPs_Bd5_A_Italia_age <- Bd5_A_Italia_merge[Bd5_A_Italia_merge[,3]%in%high_effect_SNPs_Bd5_A_Italia_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_All_A_Italia_age_full <- as.data.frame(rbind(high_effect_SNPs_Bd1_A_Italia_age, high_effect_SNPs_Bd2_A_Italia_age, high_effect_SNPs_Bd3_A_Italia_age, high_effect_SNPs_Bd4_A_Italia_age, high_effect_SNPs_Bd5_A_Italia_age))
  
  
  non_synonymous_A_Italia.vcf <- read.vcfR("./zero_fold_degenerate_SNP_Bdis326_A_Italia.vcf", verbose = FALSE )
  
  non_synonymous_Bd1_A_Italia_position <- as.numeric(non_synonymous_A_Italia.vcf@fix[non_synonymous_A_Italia.vcf@fix[,1]=="Bd1",2])
  non_synonymous_Bd1_A_Italia_age <- Bd1_A_Italia_merge[Bd1_A_Italia_merge[,3]%in%non_synonymous_Bd1_A_Italia_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd2_A_Italia_position <- as.numeric(non_synonymous_A_Italia.vcf@fix[non_synonymous_A_Italia.vcf@fix[,1]=="Bd2",2])
  non_synonymous_Bd2_A_Italia_age <- Bd2_A_Italia_merge[Bd2_A_Italia_merge[,3]%in%non_synonymous_Bd2_A_Italia_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd3_A_Italia_position <- as.numeric(non_synonymous_A_Italia.vcf@fix[non_synonymous_A_Italia.vcf@fix[,1]=="Bd3",2])
  non_synonymous_Bd3_A_Italia_age <- Bd3_A_Italia_merge[Bd3_A_Italia_merge[,3]%in%non_synonymous_Bd3_A_Italia_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd4_A_Italia_position <- as.numeric(non_synonymous_A_Italia.vcf@fix[non_synonymous_A_Italia.vcf@fix[,1]=="Bd4",2])
  non_synonymous_Bd4_A_Italia_age <- Bd4_A_Italia_merge[Bd4_A_Italia_merge[,3]%in%non_synonymous_Bd4_A_Italia_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd5_A_Italia_position <- as.numeric(non_synonymous_A_Italia.vcf@fix[non_synonymous_A_Italia.vcf@fix[,1]=="Bd5",2])
  non_synonymous_Bd5_A_Italia_age <- Bd5_A_Italia_merge[Bd5_A_Italia_merge[,3]%in%non_synonymous_Bd5_A_Italia_position, c(1:3, 15, 20)]
  
  non_synonymous_All_A_Italia_age_full <- as.data.frame(rbind(non_synonymous_Bd1_A_Italia_age, non_synonymous_Bd2_A_Italia_age, non_synonymous_Bd3_A_Italia_age, non_synonymous_Bd4_A_Italia_age, non_synonymous_Bd5_A_Italia_age))
  
}
## A East
{
  high_effect_SNPs_A_East.vcf <- read.vcfR("./Bdis326_SNPs_filtered_snpEff_high_effect.vcf", verbose = FALSE )
  
  high_effect_SNPs_Bd1_A_East_position <- as.numeric(high_effect_SNPs_A_East.vcf@fix[high_effect_SNPs_A_East.vcf@fix[,1]=="1",2])
  high_effect_SNPs_Bd1_A_East_age <- Bd1_A_East_merge[Bd1_A_East_merge[,3]%in%high_effect_SNPs_Bd1_A_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd2_A_East_position <- as.numeric(high_effect_SNPs_A_East.vcf@fix[high_effect_SNPs_A_East.vcf@fix[,1]=="2",2])
  high_effect_SNPs_Bd2_A_East_age <- Bd2_A_East_merge[Bd2_A_East_merge[,3]%in%high_effect_SNPs_Bd2_A_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd3_A_East_position <- as.numeric(high_effect_SNPs_A_East.vcf@fix[high_effect_SNPs_A_East.vcf@fix[,1]=="3",2])
  high_effect_SNPs_Bd3_A_East_age <- Bd3_A_East_merge[Bd3_A_East_merge[,3]%in%high_effect_SNPs_Bd3_A_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd4_A_East_position <- as.numeric(high_effect_SNPs_A_East.vcf@fix[high_effect_SNPs_A_East.vcf@fix[,1]=="4",2])
  high_effect_SNPs_Bd4_A_East_age <- Bd4_A_East_merge[Bd4_A_East_merge[,3]%in%high_effect_SNPs_Bd4_A_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd5_A_East_position <- as.numeric(high_effect_SNPs_A_East.vcf@fix[high_effect_SNPs_A_East.vcf@fix[,1]=="5",2])
  high_effect_SNPs_Bd5_A_East_age <- Bd5_A_East_merge[Bd5_A_East_merge[,3]%in%high_effect_SNPs_Bd5_A_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_All_A_East_age_full <- as.data.frame(rbind(high_effect_SNPs_Bd1_A_East_age, high_effect_SNPs_Bd2_A_East_age, high_effect_SNPs_Bd3_A_East_age, high_effect_SNPs_Bd4_A_East_age, high_effect_SNPs_Bd5_A_East_age))
  
  
  non_synonymous_A_East.vcf <- read.vcfR("./zero_fold_degenerate_SNP_Bdis326_A_East.vcf", verbose = FALSE )
  
  non_synonymous_Bd1_A_East_position <- as.numeric(non_synonymous_A_East.vcf@fix[non_synonymous_A_East.vcf@fix[,1]=="Bd1",2])
  non_synonymous_Bd1_A_East_age <- Bd1_A_East_merge[Bd1_A_East_merge[,3]%in%non_synonymous_Bd1_A_East_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd2_A_East_position <- as.numeric(non_synonymous_A_East.vcf@fix[non_synonymous_A_East.vcf@fix[,1]=="Bd2",2])
  non_synonymous_Bd2_A_East_age <- Bd2_A_East_merge[Bd2_A_East_merge[,3]%in%non_synonymous_Bd2_A_East_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd3_A_East_position <- as.numeric(non_synonymous_A_East.vcf@fix[non_synonymous_A_East.vcf@fix[,1]=="Bd3",2])
  non_synonymous_Bd3_A_East_age <- Bd3_A_East_merge[Bd3_A_East_merge[,3]%in%non_synonymous_Bd3_A_East_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd4_A_East_position <- as.numeric(non_synonymous_A_East.vcf@fix[non_synonymous_A_East.vcf@fix[,1]=="Bd4",2])
  non_synonymous_Bd4_A_East_age <- Bd4_A_East_merge[Bd4_A_East_merge[,3]%in%non_synonymous_Bd4_A_East_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd5_A_East_position <- as.numeric(non_synonymous_A_East.vcf@fix[non_synonymous_A_East.vcf@fix[,1]=="Bd5",2])
  non_synonymous_Bd5_A_East_age <- Bd5_A_East_merge[Bd5_A_East_merge[,3]%in%non_synonymous_Bd5_A_East_position, c(1:3, 15, 20)]
  
  non_synonymous_All_A_East_age_full <- as.data.frame(rbind(non_synonymous_Bd1_A_East_age, non_synonymous_Bd2_A_East_age, non_synonymous_Bd3_A_East_age, non_synonymous_Bd4_A_East_age, non_synonymous_Bd5_A_East_age))
  
}
## B East
{
  high_effect_SNPs_B_East.vcf <- read.vcfR("./Bdis326_SNPs_filtered_snpEff_high_effect.vcf", verbose = FALSE )
  
  high_effect_SNPs_Bd1_B_East_position <- as.numeric(high_effect_SNPs_B_East.vcf@fix[high_effect_SNPs_B_East.vcf@fix[,1]=="1",2])
  high_effect_SNPs_Bd1_B_East_age <- Bd1_B_East_merge[Bd1_B_East_merge[,3]%in%high_effect_SNPs_Bd1_B_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd2_B_East_position <- as.numeric(high_effect_SNPs_B_East.vcf@fix[high_effect_SNPs_B_East.vcf@fix[,1]=="2",2])
  high_effect_SNPs_Bd2_B_East_age <- Bd2_B_East_merge[Bd2_B_East_merge[,3]%in%high_effect_SNPs_Bd2_B_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd3_B_East_position <- as.numeric(high_effect_SNPs_B_East.vcf@fix[high_effect_SNPs_B_East.vcf@fix[,1]=="3",2])
  high_effect_SNPs_Bd3_B_East_age <- Bd3_B_East_merge[Bd3_B_East_merge[,3]%in%high_effect_SNPs_Bd3_B_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd4_B_East_position <- as.numeric(high_effect_SNPs_B_East.vcf@fix[high_effect_SNPs_B_East.vcf@fix[,1]=="4",2])
  high_effect_SNPs_Bd4_B_East_age <- Bd4_B_East_merge[Bd4_B_East_merge[,3]%in%high_effect_SNPs_Bd4_B_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_Bd5_B_East_position <- as.numeric(high_effect_SNPs_B_East.vcf@fix[high_effect_SNPs_B_East.vcf@fix[,1]=="5",2])
  high_effect_SNPs_Bd5_B_East_age <- Bd5_B_East_merge[Bd5_B_East_merge[,3]%in%high_effect_SNPs_Bd5_B_East_position, c(1:3, 15, 20)]
  
  high_effect_SNPs_All_B_East_age_full <- as.data.frame(rbind(high_effect_SNPs_Bd1_B_East_age, high_effect_SNPs_Bd2_B_East_age, high_effect_SNPs_Bd3_B_East_age, high_effect_SNPs_Bd4_B_East_age, high_effect_SNPs_Bd5_B_East_age))
  
  
  non_synonymous_B_East.vcf <- read.vcfR("./zero_fold_degenerate_SNP_Bdis326_B_East.vcf", verbose = FALSE )
  
  non_synonymous_Bd1_B_East_position <- as.numeric(non_synonymous_B_East.vcf@fix[non_synonymous_B_East.vcf@fix[,1]=="Bd1",2])
  non_synonymous_Bd1_B_East_age <- Bd1_B_East_merge[Bd1_B_East_merge[,3]%in%non_synonymous_Bd1_B_East_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd2_B_East_position <- as.numeric(non_synonymous_B_East.vcf@fix[non_synonymous_B_East.vcf@fix[,1]=="Bd2",2])
  non_synonymous_Bd2_B_East_age <- Bd2_B_East_merge[Bd2_B_East_merge[,3]%in%non_synonymous_Bd2_B_East_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd3_B_East_position <- as.numeric(non_synonymous_B_East.vcf@fix[non_synonymous_B_East.vcf@fix[,1]=="Bd3",2])
  non_synonymous_Bd3_B_East_age <- Bd3_B_East_merge[Bd3_B_East_merge[,3]%in%non_synonymous_Bd3_B_East_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd4_B_East_position <- as.numeric(non_synonymous_B_East.vcf@fix[non_synonymous_B_East.vcf@fix[,1]=="Bd4",2])
  non_synonymous_Bd4_B_East_age <- Bd4_B_East_merge[Bd4_B_East_merge[,3]%in%non_synonymous_Bd4_B_East_position, c(1:3, 15, 20)]
  
  non_synonymous_Bd5_B_East_position <- as.numeric(non_synonymous_B_East.vcf@fix[non_synonymous_B_East.vcf@fix[,1]=="Bd5",2])
  non_synonymous_Bd5_B_East_age <- Bd5_B_East_merge[Bd5_B_East_merge[,3]%in%non_synonymous_Bd5_B_East_position, c(1:3, 15, 20)]
  
  non_synonymous_All_B_East_age_full <- as.data.frame(rbind(non_synonymous_Bd1_B_East_age, non_synonymous_Bd2_B_East_age, non_synonymous_Bd3_B_East_age, non_synonymous_Bd4_B_East_age, non_synonymous_Bd5_B_East_age))
  
}


# only use 10k sites
high_effect_SNPs_All_B_West_age <- high_effect_SNPs_All_B_West_age_full[sample(1:dim(high_effect_SNPs_All_B_West_age_full)[1], size = 10000, replace = FALSE),]
non_synonymous_All_B_West_age <- non_synonymous_All_B_West_age_full[sample(1:dim(non_synonymous_All_B_West_age_full)[1], size = 10000, replace = FALSE),]
high_effect_SNPs_All_A_Italia_age <- high_effect_SNPs_All_A_Italia_age_full ## there are not more than 10k sites
non_synonymous_All_A_Italia_age <- non_synonymous_All_A_Italia_age_full[sample(1:dim(non_synonymous_All_A_Italia_age_full)[1], size = 10000, replace = FALSE),]
high_effect_SNPs_All_A_East_age <- high_effect_SNPs_All_A_East_age_full ## there are not more than 10k sites
non_synonymous_All_A_East_age <- non_synonymous_All_A_East_age_full[sample(1:dim(non_synonymous_All_A_East_age_full)[1], size = 10000, replace = FALSE),]
high_effect_SNPs_All_B_East_age <- high_effect_SNPs_All_B_East_age_full ## there are not more than 10k sites
non_synonymous_All_B_East_age <- non_synonymous_All_B_East_age_full[sample(1:dim(non_synonymous_All_B_East_age_full)[1], size = 10000, replace = FALSE),]


# downsapled neutral data
resampled_synonymous_All_age_B_West_high_effect_SNPs <- downsampling.neutral.SNP.function(synonymous_All_B_West_age, high_effect_SNPs_All_B_West_age)
resampled_synonymous_All_age_B_West_non_synonymous <- downsampling.neutral.SNP.function(synonymous_All_B_West_age, non_synonymous_All_B_West_age)
resampled_synonymous_All_age_A_Italia_high_effect_SNPs <- downsampling.neutral.SNP.function(synonymous_All_A_Italia_age, high_effect_SNPs_All_A_Italia_age)
resampled_synonymous_All_age_A_Italia_non_synonymous <- downsampling.neutral.SNP.function(synonymous_All_A_Italia_age, non_synonymous_All_A_Italia_age)
resampled_synonymous_All_age_A_East_high_effect_SNPs <- downsampling.neutral.SNP.function(synonymous_All_A_East_age, high_effect_SNPs_All_A_East_age)
resampled_synonymous_All_age_A_East_non_synonymous <- downsampling.neutral.SNP.function(synonymous_All_A_East_age, non_synonymous_All_A_East_age)
resampled_synonymous_All_age_B_East_high_effect_SNPs <- downsampling.neutral.SNP.function(synonymous_All_B_East_age, high_effect_SNPs_All_B_East_age)
resampled_synonymous_All_age_B_East_non_synonymous <- downsampling.neutral.SNP.function(synonymous_All_B_East_age, non_synonymous_All_B_East_age)


## sort data 
high_effect_SNPs_All_age_sorted_B_West <- high_effect_SNPs_All_B_West_age[order(high_effect_SNPs_All_B_West_age$PostMode), ]
resampled_synonymous_All_age_sorted_B_West_high_effect_SNPs <- resampled_synonymous_All_age_B_West_high_effect_SNPs[order(resampled_synonymous_All_age_B_West_high_effect_SNPs[,5]),]
non_synonymous_All_age_sorted_B_West <- non_synonymous_All_B_West_age[order(non_synonymous_All_B_West_age$PostMode), ]
resampled_synonymous_All_age_sorted_B_West_non_synonymous <- resampled_synonymous_All_age_B_West_non_synonymous[order(resampled_synonymous_All_age_B_West_non_synonymous[,5]),]
high_effect_SNPs_All_age_sorted_A_Italia <- high_effect_SNPs_All_A_Italia_age[order(high_effect_SNPs_All_A_Italia_age$PostMode), ]
resampled_synonymous_All_age_sorted_A_Italia_high_effect_SNPs <- resampled_synonymous_All_age_A_Italia_high_effect_SNPs[order(resampled_synonymous_All_age_A_Italia_high_effect_SNPs[,5]),]
non_synonymous_All_age_sorted_A_Italia <- non_synonymous_All_A_Italia_age[order(non_synonymous_All_A_Italia_age$PostMode), ]
resampled_synonymous_All_age_sorted_A_Italia_non_synonymous <- resampled_synonymous_All_age_A_Italia_non_synonymous[order(resampled_synonymous_All_age_A_Italia_non_synonymous[,5]),]
high_effect_SNPs_All_age_sorted_A_East <- high_effect_SNPs_All_A_East_age[order(high_effect_SNPs_All_A_East_age$PostMode), ]
resampled_synonymous_All_age_sorted_A_East_high_effect_SNPs <- resampled_synonymous_All_age_A_East_high_effect_SNPs[order(resampled_synonymous_All_age_A_East_high_effect_SNPs[,5]),]
non_synonymous_All_age_sorted_A_East <- non_synonymous_All_A_East_age[order(non_synonymous_All_A_East_age$PostMode), ]
resampled_synonymous_All_age_sorted_A_East_non_synonymous <- resampled_synonymous_All_age_A_East_non_synonymous[order(resampled_synonymous_All_age_A_East_non_synonymous[,5]),]
high_effect_SNPs_All_age_sorted_B_East <- high_effect_SNPs_All_B_East_age[order(high_effect_SNPs_All_B_East_age$PostMode), ]
resampled_synonymous_All_age_sorted_B_East_high_effect_SNPs <- resampled_synonymous_All_age_B_East_high_effect_SNPs[order(resampled_synonymous_All_age_B_East_high_effect_SNPs[,5]),]
non_synonymous_All_age_sorted_B_East <- non_synonymous_All_B_East_age[order(non_synonymous_All_B_East_age$PostMode), ]
resampled_synonymous_All_age_sorted_B_East_non_synonymous <- resampled_synonymous_All_age_B_East_non_synonymous[order(resampled_synonymous_All_age_B_East_non_synonymous[,5]),]


## get delta frequency:
delta_freq_decile_All_resample_B_West_high_effect_SNPs <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_B_West_high_effect_SNPs, high_effect_SNPs_All_age_sorted_B_West, my_bonderies_B_West)
delta_freq_decile_All_resample_B_West_non_synonymous <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_B_West_non_synonymous, non_synonymous_All_age_sorted_B_West, my_bonderies_B_West)
delta_freq_decile_All_resample_A_Italia_high_effect_SNPs <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_A_Italia_high_effect_SNPs, high_effect_SNPs_All_age_sorted_A_Italia, my_bonderies_A_Italia)
delta_freq_decile_All_resample_A_Italia_non_synonymous <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_A_Italia_non_synonymous, non_synonymous_All_age_sorted_A_Italia, my_bonderies_A_Italia)
delta_freq_decile_All_resample_A_East_high_effect_SNPs <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_A_East_high_effect_SNPs, high_effect_SNPs_All_age_sorted_A_East, my_bonderies_A_East)
delta_freq_decile_All_resample_A_East_non_synonymous <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_A_East_non_synonymous, non_synonymous_All_age_sorted_A_East, my_bonderies_A_East)
delta_freq_decile_All_resample_B_East_high_effect_SNPs <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_B_East_high_effect_SNPs, high_effect_SNPs_All_age_sorted_B_East, my_bonderies_B_East)
delta_freq_decile_All_resample_B_East_non_synonymous <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_B_East_non_synonymous, non_synonymous_All_age_sorted_B_East, my_bonderies_B_East)


## bootstrapping
delta_freq_decile_All_boot_B_West_high_effect_SNPs <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- high_effect_SNPs_All_B_West_age[sample(1:dim(high_effect_SNPs_All_B_West_age)[1], dim(high_effect_SNPs_All_B_West_age)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(synonymous_All_B_West_age, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_West)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_B_West_non_synonymous <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- non_synonymous_All_age_sorted_B_West[sample(1:dim(non_synonymous_All_age_sorted_B_West)[1], dim(non_synonymous_All_age_sorted_B_West)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(synonymous_All_B_West_age, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_West)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_A_Italia_high_effect_SNPs <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- high_effect_SNPs_All_A_Italia_age[sample(1:dim(high_effect_SNPs_All_A_Italia_age)[1], dim(high_effect_SNPs_All_A_Italia_age)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(synonymous_All_A_Italia_age, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_Italia)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_A_Italia_non_synonymous <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- non_synonymous_All_age_sorted_A_Italia[sample(1:dim(non_synonymous_All_age_sorted_A_Italia)[1], dim(non_synonymous_All_age_sorted_A_Italia)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(synonymous_All_A_Italia_age, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_Italia)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_A_East_high_effect_SNPs <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- high_effect_SNPs_All_A_East_age[sample(1:dim(high_effect_SNPs_All_A_East_age)[1], dim(high_effect_SNPs_All_A_East_age)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(synonymous_All_A_East_age, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_East)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_A_East_non_synonymous <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- non_synonymous_All_age_sorted_A_East[sample(1:dim(non_synonymous_All_age_sorted_A_East)[1], dim(non_synonymous_All_age_sorted_A_East)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(synonymous_All_A_East_age, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_East)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_B_East_high_effect_SNPs <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- high_effect_SNPs_All_B_East_age[sample(1:dim(high_effect_SNPs_All_B_East_age)[1], dim(high_effect_SNPs_All_B_East_age)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(synonymous_All_B_East_age, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_East)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_B_East_non_synonymous <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- non_synonymous_All_age_sorted_B_East[sample(1:dim(non_synonymous_All_age_sorted_B_East)[1], dim(non_synonymous_All_age_sorted_B_East)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(synonymous_All_B_East_age, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_East)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)












