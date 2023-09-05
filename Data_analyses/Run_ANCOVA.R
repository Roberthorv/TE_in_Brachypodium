## Run ANOVA based on the results of the iHS analyses

## read in iHS regions
high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab <- read.table("./high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab.bed", sep = "\t", header = FALSE)

##select the right TEs frequencies based on clades
TEs_in_high_iHS_regions_all_clades_pos_sel_less_5_perc <- read.table("./TEs_in_high_iHS_regions_all_clades_pos_sel_less_5_perc.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

## make the table
Fixed_TE_proportion_in_iHS <- matrix(NA, nrow = length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))*3*4*dim(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab)[1] , ncol = 7)
colnames(Fixed_TE_proportion_in_iHS) <- c("Region", "Proportion_of_fixed_TEs", "Total_number_of_TEs", "Hihg_iHS", "TE_Superfamily", "TE_age", "Clade")

## get unique region names
for (i in 1:dim(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab)[1]) {
  Fixed_TE_proportion_in_iHS[1:(length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))*3*4)+(i-1)*(length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))*3*4),1] <- paste(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[i,1:3], collapse = "_")
}

## add clades
Fixed_TE_proportion_in_iHS[,7] <- c(rep("A_East", length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))*3), 
                                    rep("A_Italy", length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))*3), 
                                    rep("B_East", length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))*3), 
                                    rep("B_West", length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))*3))

## add age
Fixed_TE_proportion_in_iHS[,6] <- c(rep("Young", length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))), 
                                    rep("Intermediate", length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))), 
                                    rep("Old", length(unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7]))))

## add TE superfamily
Fixed_TE_proportion_in_iHS[,5] <- unique(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,7])

## add iHS regions
Fixed_TE_proportion_in_iHS[,4] <- 0
for (i in 1:dim(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab)[1]) {
  my_region <- paste(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[i,1:3], collapse = "_") 
  my_clade <- unlist(strsplit(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[i,4], split = ","))
  Fixed_TE_proportion_in_iHS[Fixed_TE_proportion_in_iHS[,1] == my_region & Fixed_TE_proportion_in_iHS[,7] %in% my_clade,4] <- 1
}

## add TE proportion and total number
## read in clade data 
clade_table <- read.table("./bdis_333_clade.csv", sep = ",", header = TRUE)
names <- gsub(pattern = "\\.", replacement = "-", gsub(pattern = "X", replacement = "", names(Full_clean_TIP_TAP_Sample_table_326_samples[11:336])))

## get names of samples in each clade:
A_East_clade_names <- c()
for (n in names) {
  if (clade_table[clade_table$names %in% n | clade_table$altern_names %in% n,11] == "A_East") {
    A_East_clade_names <- c(A_East_clade_names, n)
  }
}

A_Italia_clade_names <- c()
for (n in names) {
  if (clade_table[clade_table$names %in% n | clade_table$altern_names %in% n,11] == "A_Italia") {
    A_Italia_clade_names <- c(A_Italia_clade_names, n)
  }
}

B_West_clade_names <- c()
for (n in names) {
  if (clade_table[clade_table$names %in% n | clade_table$altern_names %in% n,11] == "B_West") {
    B_West_clade_names <- c(B_West_clade_names, n)
  }
}

B_East_clade_names <- c()
for (n in names) {
  if (clade_table[clade_table$names %in% n | clade_table$altern_names %in% n,11] == "B_East") {
    B_East_clade_names <- c(B_East_clade_names, n)
  }
}

C_clade_names <- c()
for (n in names) {
  if (clade_table[clade_table$names %in% n | clade_table$altern_names %in% n,11] == "C") {
    C_clade_names <- c(C_clade_names, n)
  }
}

## find multialleic TEs:
multiallelic_TE_ID <- Full_clean_TIP_TAP_Sample_table_326_samples[grep(pattern = ",", Full_clean_TIP_TAP_Sample_table_326_samples[,6]),4]

## remove multialleleic TEs
Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples <- Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,]

## All clade age
{
  Bd1_All_clusters_Clock_J_age_estimate <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd1_full_A_B_cluster_Clock_J_age_estimate.sites2.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE)
  colnames(Bd1_All_clusters_Clock_J_age_estimate) <- c("MarkerID", "Clock", "N_Concordant", "N_Discordant", "PostMode")
  Bd1_All_clusters.marker <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd1_full_A_B_cluster.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd2_All_clusters_Clock_J_age_estimate <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd2_full_A_B_cluster_Clock_J_age_estimate.sites2.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE)
  colnames(Bd2_All_clusters_Clock_J_age_estimate) <- c("MarkerID", "Clock", "N_Concordant", "N_Discordant", "PostMode")
  Bd2_All_clusters.marker <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd2_full_A_B_cluster.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd3_All_clusters_Clock_J_age_estimate <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd3_full_A_B_cluster_Clock_J_age_estimate.sites2.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE)
  colnames(Bd3_All_clusters_Clock_J_age_estimate) <- c("MarkerID", "Clock", "N_Concordant", "N_Discordant", "PostMode")
  Bd3_All_clusters.marker <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd3_full_A_B_cluster.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd4_All_clusters_Clock_J_age_estimate <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd4_full_A_B_cluster_Clock_J_age_estimate.sites2.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE)
  colnames(Bd4_All_clusters_Clock_J_age_estimate) <- c("MarkerID", "Clock", "N_Concordant", "N_Discordant", "PostMode")
  Bd4_All_clusters.marker <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd4_full_A_B_cluster.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd5_All_clusters_Clock_J_age_estimate <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd5_full_A_B_cluster_Clock_J_age_estimate.sites2.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE)
  colnames(Bd5_All_clusters_Clock_J_age_estimate) <- c("MarkerID", "Clock", "N_Concordant", "N_Discordant", "PostMode")
  Bd5_All_clusters.marker <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/GEVA/All_clusters/Bd5_full_A_B_cluster.marker.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
  
  Bd1_All_clusters_merge <- merge(Bd1_All_clusters.marker, Bd1_All_clusters_Clock_J_age_estimate, by = "MarkerID")
  Bd2_All_clusters_merge <- merge(Bd2_All_clusters.marker, Bd2_All_clusters_Clock_J_age_estimate, by = "MarkerID")
  Bd3_All_clusters_merge <- merge(Bd3_All_clusters.marker, Bd3_All_clusters_Clock_J_age_estimate, by = "MarkerID")
  Bd4_All_clusters_merge <- merge(Bd4_All_clusters.marker, Bd4_All_clusters_Clock_J_age_estimate, by = "MarkerID")
  Bd5_All_clusters_merge <- merge(Bd5_All_clusters.marker, Bd5_All_clusters_Clock_J_age_estimate, by = "MarkerID")
  
  All_All_clusters_age_estimate <- rbind(Bd1_All_clusters_merge, Bd2_All_clusters_merge, Bd3_All_clusters_merge, Bd4_All_clusters_merge, Bd5_All_clusters_merge)
  
}


## get relevent TEs, the freqeucies and ages
## per clades
TEs_in_regions_of_interest_A_East_ID <- c()
for (t in 1:dim(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab)[1]) {
  my_TE_IDs <- Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,1] == high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,1] &
                                                                      Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,2] <= high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,3] &
                                                                       Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,3] >= high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,2], 4]
  TEs_in_regions_of_interest_A_East_ID <- c(TEs_in_regions_of_interest_A_East_ID, my_TE_IDs)
}

TEs_in_regions_of_interest_A_East_raw <- matrix(NA, nrow = length(TEs_in_regions_of_interest_A_East_ID), ncol = 7)
TEs_in_regions_of_interest_A_East_raw[,4] <- TEs_in_regions_of_interest_A_East_ID

for (t in TEs_in_regions_of_interest_A_East_ID) {
  TEs_in_regions_of_interest_A_East_raw[TEs_in_regions_of_interest_A_East_raw[,4] == t, 1:3] <- unlist(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 1:3])
  TEs_in_regions_of_interest_A_East_raw[TEs_in_regions_of_interest_A_East_raw[,4] == t, 5] <- unlist(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 7])
  DA_countr <- sum(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples) %in% A_East_clade_names], na.rm = TRUE)
  total_A_count <- sum(!is.na(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples) %in% A_East_clade_names]))
  TEs_in_regions_of_interest_A_East_raw[TEs_in_regions_of_interest_A_East_raw[,4] == t, 6] <- DA_countr/total_A_count
  my_age <- All_All_clusters_age_estimate[All_All_clusters_age_estimate$Label == t, 18]
  if (length(my_age) == 1) {
    TEs_in_regions_of_interest_A_East_raw[TEs_in_regions_of_interest_A_East_raw[,4] == t, 7] <-  my_age
  }

}
colnames(TEs_in_regions_of_interest_A_East_raw) <- c("CHR", "Start", "End", "ID", "TE_Superfamily", "Freq", "Age")
TEs_in_regions_of_interest_A_East <- TEs_in_regions_of_interest_A_East_raw[TEs_in_regions_of_interest_A_East_raw[,6] != "0" & !is.na(TEs_in_regions_of_interest_A_East_raw[,7]),]


TEs_in_regions_of_interest_A_Italia_ID <- c()
for (t in 1:dim(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab)[1]) {
  my_TE_IDs <- Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,1] == high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,1] &
                                                                       Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,2] <= high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,3] &
                                                                       Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,3] >= high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,2], 4]
  TEs_in_regions_of_interest_A_Italia_ID <- c(TEs_in_regions_of_interest_A_Italia_ID, my_TE_IDs)
}

TEs_in_regions_of_interest_A_Italia_raw <- matrix(NA, nrow = length(TEs_in_regions_of_interest_A_Italia_ID), ncol = 7)
TEs_in_regions_of_interest_A_Italia_raw[,4] <- TEs_in_regions_of_interest_A_Italia_ID

for (t in TEs_in_regions_of_interest_A_Italia_ID) {
  TEs_in_regions_of_interest_A_Italia_raw[TEs_in_regions_of_interest_A_Italia_raw[,4] == t, 1:3] <- unlist(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 1:3])
  TEs_in_regions_of_interest_A_Italia_raw[TEs_in_regions_of_interest_A_Italia_raw[,4] == t, 5] <- unlist(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 7])
  DA_countr <- sum(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples) %in% A_Italia_clade_names], na.rm = TRUE)
  total_A_count <- sum(!is.na(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples) %in% A_Italia_clade_names]))
  TEs_in_regions_of_interest_A_Italia_raw[TEs_in_regions_of_interest_A_Italia_raw[,4] == t, 6] <- DA_countr/total_A_count
  my_age <- All_All_clusters_age_estimate[All_All_clusters_age_estimate$Label == t, 18]
  if (length(my_age) == 1) {
    TEs_in_regions_of_interest_A_Italia_raw[TEs_in_regions_of_interest_A_Italia_raw[,4] == t, 7] <-  my_age
  }
  
}
colnames(TEs_in_regions_of_interest_A_Italia_raw) <- c("CHR", "Start", "End", "ID", "TE_Superfamily", "Freq", "Age")
TEs_in_regions_of_interest_A_Italia <- TEs_in_regions_of_interest_A_Italia_raw[TEs_in_regions_of_interest_A_Italia_raw[,6] != "0" & !is.na(TEs_in_regions_of_interest_A_Italia_raw[,7]),]


TEs_in_regions_of_interest_B_West_ID <- c()
for (t in 1:dim(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab)[1]) {
  my_TE_IDs <- Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,1] == high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,1] &
                                                                       Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,2] <= high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,3] &
                                                                       Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,3] >= high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,2], 4]
  TEs_in_regions_of_interest_B_West_ID <- c(TEs_in_regions_of_interest_B_West_ID, my_TE_IDs)
}

TEs_in_regions_of_interest_B_West_raw <- matrix(NA, nrow = length(TEs_in_regions_of_interest_B_West_ID), ncol = 7)
TEs_in_regions_of_interest_B_West_raw[,4] <- TEs_in_regions_of_interest_B_West_ID

for (t in TEs_in_regions_of_interest_B_West_ID) {
  TEs_in_regions_of_interest_B_West_raw[TEs_in_regions_of_interest_B_West_raw[,4] == t, 1:3] <- unlist(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 1:3])
  TEs_in_regions_of_interest_B_West_raw[TEs_in_regions_of_interest_B_West_raw[,4] == t, 5] <- unlist(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 7])
  DA_countr <- sum(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples) %in% B_West_clade_names], na.rm = TRUE)
  total_A_count <- sum(!is.na(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples) %in% B_West_clade_names]))
  TEs_in_regions_of_interest_B_West_raw[TEs_in_regions_of_interest_B_West_raw[,4] == t, 6] <- DA_countr/total_A_count
  my_age <- All_All_clusters_age_estimate[All_All_clusters_age_estimate$Label == t, 18]
  if (length(my_age) == 1) {
    TEs_in_regions_of_interest_B_West_raw[TEs_in_regions_of_interest_B_West_raw[,4] == t, 7] <-  my_age
  }
  
}
colnames(TEs_in_regions_of_interest_B_West_raw) <- c("CHR", "Start", "End", "ID", "TE_Superfamily", "Freq", "Age")
TEs_in_regions_of_interest_B_West <- TEs_in_regions_of_interest_B_West_raw[TEs_in_regions_of_interest_B_West_raw[,6] != "0" & !is.na(TEs_in_regions_of_interest_B_West_raw[,7]),]


TEs_in_regions_of_interest_B_East_ID <- c()
for (t in 1:dim(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab)[1]) {
  my_TE_IDs <- Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,1] == high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,1] &
                                                                       Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,2] <= high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,3] &
                                                                       Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,3] >= high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[t,2], 4]
  TEs_in_regions_of_interest_B_East_ID <- c(TEs_in_regions_of_interest_B_East_ID, my_TE_IDs)
}

TEs_in_regions_of_interest_B_East_raw <- matrix(NA, nrow = length(TEs_in_regions_of_interest_B_East_ID), ncol = 7)
TEs_in_regions_of_interest_B_East_raw[,4] <- TEs_in_regions_of_interest_B_East_ID

for (t in TEs_in_regions_of_interest_B_East_ID) {
  TEs_in_regions_of_interest_B_East_raw[TEs_in_regions_of_interest_B_East_raw[,4] == t, 1:3] <- unlist(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 1:3])
  TEs_in_regions_of_interest_B_East_raw[TEs_in_regions_of_interest_B_East_raw[,4] == t, 5] <- unlist(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 7])
  DA_countr <- sum(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples) %in% B_East_clade_names], na.rm = TRUE)
  total_A_count <- sum(!is.na(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Biallelic_Full_clean_TIP_TAP_Sample_table_326_samples) %in% B_East_clade_names]))
  TEs_in_regions_of_interest_B_East_raw[TEs_in_regions_of_interest_B_East_raw[,4] == t, 6] <- DA_countr/total_A_count
  my_age <- All_All_clusters_age_estimate[All_All_clusters_age_estimate$Label == t, 18]
  if (length(my_age) == 1) {
    TEs_in_regions_of_interest_B_East_raw[TEs_in_regions_of_interest_B_East_raw[,4] == t, 7] <-  my_age
  }
  
}
colnames(TEs_in_regions_of_interest_B_East_raw) <- c("CHR", "Start", "End", "ID", "TE_Superfamily", "Freq", "Age")
TEs_in_regions_of_interest_B_East <- TEs_in_regions_of_interest_B_East_raw[TEs_in_regions_of_interest_B_East_raw[,6] != "0" & !is.na(TEs_in_regions_of_interest_B_East_raw[,7]),]




## age splits: 
## AB split: 52000 generation ago
## B split: 22000 generation ago
## A split: 18000 generation ago


## Old > 60000 generation 
## Intermediate < 60000 & > 10000 generation 
## Young < 10000 generation 

for (i in 1:dim(Fixed_TE_proportion_in_iHS)[1]) {
  my_TE_Superfamily <- Fixed_TE_proportion_in_iHS[i,5]
  my_age <- Fixed_TE_proportion_in_iHS[i,6]
  my_clade <- Fixed_TE_proportion_in_iHS[i,7]
  if(my_clade == "A_Italy"){ my_clade <- "A_Italia"}
  my_scaff <- unlist(strsplit(Fixed_TE_proportion_in_iHS[i,1], split = "_"))[1]
  my_start <- unlist(strsplit(Fixed_TE_proportion_in_iHS[i,1], split = "_"))[2]
  my_end <- unlist(strsplit(Fixed_TE_proportion_in_iHS[i,1], split = "_"))[3]
  my_tab_name <- paste0("TEs_in_regions_of_interest_", my_clade)
  
  
  my_run_tab <- get(my_tab_name)[get(my_tab_name)[,1] == my_scaff &
                                   as.numeric(get(my_tab_name)[,2]) >= as.numeric(my_start) &
                                   as.numeric(get(my_tab_name)[,3]) <= as.numeric(my_end) &
                                   get(my_tab_name)[,5] == my_TE_Superfamily,]
  
  if (class(my_run_tab) %in% c("matrix", "array" ) ) {
    if (my_age == "Young") {
      my_run_tab_age <- my_run_tab[my_run_tab[, 7] < 10000,]
    } else if (my_age == "Intermediate") {
      my_run_tab_age <- my_run_tab[my_run_tab[, 7] >= 10000 & my_run_tab[, 7] < 60000,]
    } else {
      my_run_tab_age <- my_run_tab[my_run_tab[, 7] >= 60000,]
    }
    
    if (class(my_run_tab_age) %in% c("matrix", "array" ) ) {
      Fixed_TE_proportion_in_iHS[i,3] <- dim(my_run_tab_age)[1]
      Fixed_TE_proportion_in_iHS[i,2] <- sum(my_run_tab_age[,6] == 1 )/dim(my_run_tab_age)[1]
    } else {
      if (length(my_run_tab_age) == 0) {
        Fixed_TE_proportion_in_iHS[i,3] <- 0
        Fixed_TE_proportion_in_iHS[i,2] <- "NaN"
      } else {
        Fixed_TE_proportion_in_iHS[i,3] <- 1
        Fixed_TE_proportion_in_iHS[i,2] <- sum(my_run_tab_age[6] == 1 )/1
      }
    }
    
  } else {
    if (my_age == "Young") {
      my_run_tab_age <- my_run_tab[my_run_tab[7] < 10000]
    } else if (my_age == "Intermediate") {
      my_run_tab_age <- my_run_tab[my_run_tab[7] >= 10000 & my_run_tab[7] < 60000]
    } else {
      my_run_tab_age <- my_run_tab[my_run_tab[7] >= 60000]
    }
    
    if (length(my_run_tab_age) == 0) {
      Fixed_TE_proportion_in_iHS[i,3] <- 0
      Fixed_TE_proportion_in_iHS[i,2] <- "NaN"
    } else {
      Fixed_TE_proportion_in_iHS[i,3] <- 1
      Fixed_TE_proportion_in_iHS[i,2] <- sum(my_run_tab_age[6] == 1 )/1
    }

  }
  
 
  print(i/dim(Fixed_TE_proportion_in_iHS)[1])
  
}


## write.table(Fixed_TE_proportion_in_iHS, file = "/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/R_Input/Fixed_TE_proportion_in_iHS.txt", quote = FALSE, sep = "\t")

## read in data table
Fixed_TE_proportion_in_iHS_df <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/R_Input/Fixed_TE_proportion_in_iHS.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

Fixed_TE_proportion_in_iHS_df[1:5,]
Fixed_TE_proportion_in_iHS_df_b <- Fixed_TE_proportion_in_iHS_df
Fixed_TE_proportion_in_iHS_df_b$Number_of_fixed_TEs <- Fixed_TE_proportion_in_iHS_df$Proportion_of_fixed_TEs * Fixed_TE_proportion_in_iHS_df$Total_number_of_TEs
Fixed_TE_proportion_in_iHS_df_b[1:5,]
Fixed_TE_proportion_in_iHS_df_b[is.na(Fixed_TE_proportion_in_iHS_df_b$Number_of_fixed_TEs),]$Number_of_fixed_TEs <- 0
Fixed_TE_proportion_in_iHS_df_b$Hihg_iHS <- as.factor(Fixed_TE_proportion_in_iHS_df_b$Hihg_iHS)

## run an ANCOVA:
library(car)
ancova_model_fix_number <- aov(Number_of_fixed_TEs ~ Total_number_of_TEs + Hihg_iHS + Clade + TE_Superfamily + Region + TE_age, data = Fixed_TE_proportion_in_iHS_df_b)
summary(ancova_model_fix_number)
Anova(ancova_model_fix_number, type="III")


## the residuals of the model are ruffly normally distributed
hist(resid(aov(Number_of_fixed_TEs ~ Total_number_of_TEs + Hihg_iHS + Clade + TE_Superfamily + Region + TE_age, data = Fixed_TE_proportion_in_iHS_df_b)), breaks = seq(-20,15,0.5))



## Hihg iHS is not effecting the number of fixed TEs
## But the tree most important factors based on the Sum sq are:
## Total_number_of_TEs, TE_Superfamily and Clade

cor(Fixed_TE_proportion_in_iHS_df_b$Total_number_of_TEs, Fixed_TE_proportion_in_iHS_df_b$Number_of_fixed_TEs)
## the more TEs there are the more are fixed

library(multcomp)
postHocs_fix_number_TE_Superfamily <- glht(ancova_model_fix_number, linfct = mcp(TE_Superfamily = "Tukey"))
summary(postHocs_fix_number_TE_Superfamily)
## TEs are behaving differently in different families

postHocs_fix_number_Clade <- glht(ancova_model_fix_number, linfct = mcp(Clade = "Tukey"))
summary(postHocs_fix_number_Clade)
## the number of fixed TEs correlates with Ne

postHocs_fix_number_iHS <- glht(ancova_model_fix_number, linfct = mcp(Hihg_iHS = "Tukey"))
summary(postHocs_fix_number_iHS)





##################################################################################
##############           Analyes the freqeuncy of the TEs           ##############
##################################################################################


## make a table with each TEs per clade per line

## make the table
TE_ID_LIST <- names(table(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,4])[table(Full_clean_TIP_TAP_Sample_table_326_samples[!Full_clean_TIP_TAP_Sample_table_326_samples[,4] %in% multiallelic_TE_ID,4]) == 1])
Frequency_TE_in_iHS <- matrix(NA, nrow = length(TE_ID_LIST)*4 , ncol = 7)
colnames(Frequency_TE_in_iHS) <- c("TE_ID", "Frequency", "iHS", "TE_Superfamily", "Age", "Clade", "Regoin")

## add clades
Frequency_TE_in_iHS[,6] <- c(rep("A_East", length(TE_ID_LIST)), 
                                    rep("A_Italy", length(TE_ID_LIST)), 
                                    rep("B_East", length(TE_ID_LIST)), 
                                    rep("B_West", length(TE_ID_LIST)))

## add TE ID
Frequency_TE_in_iHS[,1] <- c(rep(TE_ID_LIST, 4))

## add rest
for (t in TE_ID_LIST) {
  for (c in c("A_East", "A_Italy", "B_East", "B_West")) {
    Frequency_TE_in_iHS[Frequency_TE_in_iHS[,1] == t & Frequency_TE_in_iHS[,6] == c, 4] <- Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 7]
    if(t %in% All_All_clusters_age_estimate$Label){
      Frequency_TE_in_iHS[Frequency_TE_in_iHS[,1] == t & Frequency_TE_in_iHS[,6] == c, 5] <- All_All_clusters_age_estimate[All_All_clusters_age_estimate$Label == t, 18]
    } else {Frequency_TE_in_iHS[Frequency_TE_in_iHS[,1] == t & Frequency_TE_in_iHS[,6] == c, 5] <- NA}
    my_Scaf <- Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 1]
    my_start <- Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 2]
    my_end <- Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, 3]
    my_iHS_reion_clades <- unlist(strsplit(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[
      high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[,1] == my_Scaf &
      high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[,2] <=  my_end &
      high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[,3] >= my_start,4], split = ","))
    if (c %in% my_iHS_reion_clades) {
      Frequency_TE_in_iHS[Frequency_TE_in_iHS[,1] == t & Frequency_TE_in_iHS[,6] == c, 3] <- "High"
    } else {Frequency_TE_in_iHS[Frequency_TE_in_iHS[,1] == t & Frequency_TE_in_iHS[,6] == c, 3] <- "Average"}
    Frequency_TE_in_iHS[Frequency_TE_in_iHS[,1] == t & Frequency_TE_in_iHS[,6] == c, 7] <- paste(high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[
      high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[,1] == my_Scaf &
        high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[,2] <=  my_end &
        high_iHS_regions_all_clades_pos_sel_less_5_perc_ID_tab[,3] >= my_start, 1:3], collapse = "_")
    my_data_name <- gsub(pattern = "Italy", replacement = "Italia", paste(c, "_clade_names", sep = ""))
    DA_countr <- sum(Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Full_clean_TIP_TAP_Sample_table_326_samples) %in% get(my_data_name)], na.rm = TRUE)
    total_A_count <- sum(!is.na(Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples[,4] == t, names(Full_clean_TIP_TAP_Sample_table_326_samples) %in% get(my_data_name)]))
    Frequency_TE_in_iHS[Frequency_TE_in_iHS[,1] == t & Frequency_TE_in_iHS[,6] == c, 2] <- DA_countr/total_A_count
  }
  print(sum(1:length(TE_ID_LIST) * c(TE_ID_LIST == t))/length(TE_ID_LIST))
}

## select only TEs in high iHS regions:
Frequency_TE_in_high_iHS <- Frequency_TE_in_iHS[grepl(pattern = "Bd", Frequency_TE_in_iHS[,7]),]

## write.table(Frequency_TE_in_high_iHS, file = "/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/R_Input/Frequency_TE_in_high_iHS.txt", quote = FALSE, sep = "\t")

## read in data table
Frequency_TE_in_high_iHS_df <- read.table("/Users/roberthorvath/Desktop/Projects/1) TEs_in_Brachypodium/4) Results/R_Input/Frequency_TE_in_high_iHS.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

## select only polymorphic TEs
Frequency_poly_TE_in_high_iHS_df <- Frequency_TE_in_high_iHS_df[Frequency_TE_in_high_iHS_df[,2] != 0 & Frequency_TE_in_high_iHS_df[,2] != 1,]

## bin age as was done previously:
Frequency_poly_TE_in_high_iHS_age_bin_df <- Frequency_poly_TE_in_high_iHS_df

## Old > 60000 generation 
## Intermediate < 60000 & > 10000 generation 
## Young < 10000 generation 

for (a in 1:dim(Frequency_poly_TE_in_high_iHS_age_bin_df)[1]) {
  if (is.na(Frequency_poly_TE_in_high_iHS_df[a,]$Age)) {
    next
  } else if (Frequency_poly_TE_in_high_iHS_df[a,]$Age < 10000) {
    Frequency_poly_TE_in_high_iHS_age_bin_df[a,]$Age <- "Young"
  } else if (Frequency_poly_TE_in_high_iHS_df[a,]$Age < 60000) {
    Frequency_poly_TE_in_high_iHS_age_bin_df[a,]$Age <- "Intermediate"
  } else {
    Frequency_poly_TE_in_high_iHS_age_bin_df[a,]$Age <- "Old"
  }
}

Frequency_poly_TE_in_high_iHS_age_bin_df$Age <- as.factor(Frequency_poly_TE_in_high_iHS_age_bin_df$Age)

## run an ANCOVA:
library(car)
ancova_model_poly_TE_age_bin <- aov(Frequency ~ iHS + Clade + TE_Superfamily + Regoin + Age, data = Frequency_poly_TE_in_high_iHS_age_bin_df)
summary(ancova_model_poly_TE_age_bin)
Anova(ancova_model_poly_TE_age_bin, type="III")

## the residuals of the model are ruffly normally distributed
hist(resid(aov(Frequency ~ iHS + Clade + TE_Superfamily + Regoin + Age, data = Frequency_poly_TE_in_high_iHS_age_bin_df)), breaks = seq(-1.5,1.5,0.1))


## Hihg iHS is not effecting ...
## But the tree most important factors based on the Sum sq are:
## TE_Superfamily, Clade and Region

library(multcomp)
postHocs_poly_TE_age_bin_TE_Superfamily <- glht(ancova_model_poly_TE_age_bin, linfct = mcp(TE_Superfamily = "Tukey"))
summary(postHocs_poly_TE_age_bin_TE_Superfamily)
## TEs are behaving differently in different families

postHocs_poly_TE_age_bin_Clade <- glht(ancova_model_poly_TE_age_bin, linfct = mcp(Clade = "Tukey"))
summary(postHocs_poly_TE_age_bin_Clade)
## the patter here is not as clear as before with the fixed TEs

postHocs_poly_TE_age_bin_Age <- glht(ancova_model_poly_TE_age_bin, linfct = mcp(Age = "Tukey"))
summary(postHocs_poly_TE_age_bin_Age)
## only young TEs are less frequent than intermediate TEs

postHocs_poly_TE_age_bin_Region <- glht(ancova_model_poly_TE_age_bin, linfct = mcp(Regoin = "Tukey"))
summary(postHocs_poly_TE_age_bin_Region)
## to many variables...

postHocs_poly_TE_age_bin_iHS <- glht(ancova_model_poly_TE_age_bin, linfct = mcp(iHS = "Tukey"))
summary(postHocs_poly_TE_age_bin_iHS)
## 



