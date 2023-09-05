## split data
TEs_in_gene_and_1kB_surroundings_B_West <- read.table("./TEs_in_gene_and_1kB_surroundings_B_West", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
TEs_in_gene_and_1kB_surroundings_A_Italia <- read.table("./TEs_in_gene_and_1kB_surroundings_A_Italia", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
TEs_in_gene_and_1kB_surroundings_A_East <- read.table("./TEs_in_gene_and_1kB_surroundings_A_East", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
TEs_in_gene_and_1kB_surroundings_B_East <- read.table("./TEs_in_gene_and_1kB_surroundings_B_East", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

TEs_in_gene_and_5kB_surroundings_B_West <- read.table("./TEs_in_gene_and_5kB_surroundings_B_West", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
TEs_in_gene_and_5kB_surroundings_A_Italia <- read.table("./TEs_in_gene_and_5kB_surroundings_A_Italia", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
TEs_in_gene_and_5kB_surroundings_A_East <- read.table("./TEs_in_gene_and_5kB_surroundings_A_East", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
TEs_in_gene_and_5kB_surroundings_B_East <- read.table("./TEs_in_gene_and_5kB_surroundings_B_East", sep = "\t", header = FALSE, stringsAsFactors = FALSE)


## get the different list of TEs:
Retrotransposons_All_age_B_East_in_gene_and_1kB <- All_B_East_merge[All_B_East_merge$Label %in% Retrotransposon_ID & !All_B_East_merge$Position %in% not_derived_TE_All_B_East[,2] & All_B_East_merge$Label %in% TEs_in_gene_and_1kB_surroundings_B_East[,4],c(1:4,15,20)]
Retrotransposons_All_age_A_Italia_in_gene_and_1kB <- All_A_Italia_merge[All_A_Italia_merge$Label %in% Retrotransposon_ID & !All_A_Italia_merge$Position %in% not_derived_TE_All_A_Italia[,2] & All_A_Italia_merge$Label %in% TEs_in_gene_and_1kB_surroundings_A_Italia[,4],c(1:4,15,20)]
Retrotransposons_All_age_A_East_in_gene_and_1kB <- All_A_East_merge[All_A_East_merge$Label %in% Retrotransposon_ID & !All_A_East_merge$Position %in% not_derived_TE_All_A_East[,2] & All_A_East_merge$Label %in% TEs_in_gene_and_1kB_surroundings_A_East[,4],c(1:4,15,20)]
Retrotransposons_All_age_B_West_in_gene_and_1kB <- All_B_West_merge[All_B_West_merge$Label %in% Retrotransposon_ID & !All_B_West_merge$Position %in% not_derived_TE_All_B_West[,2] & All_B_West_merge$Label %in% TEs_in_gene_and_1kB_surroundings_B_West[,4],c(1:4,15,20)]

Retrotransposons_All_age_B_East_1kB_and_5kB <- All_B_East_merge[All_B_East_merge$Label %in% Retrotransposon_ID & !All_B_East_merge$Position %in% not_derived_TE_All_B_East[,2] & !All_B_East_merge$Label %in% TEs_in_gene_and_1kB_surroundings_B_East[,4] & All_B_East_merge$Label %in% TEs_in_gene_and_5kB_surroundings_B_East[,4],c(1:4,15,20)]
Retrotransposons_All_age_A_Italia_1kB_and_5kB <- All_A_Italia_merge[All_A_Italia_merge$Label %in% Retrotransposon_ID & !All_A_Italia_merge$Position %in% not_derived_TE_All_A_Italia[,2] & !All_A_Italia_merge$Label %in% TEs_in_gene_and_1kB_surroundings_A_Italia[,4] & All_A_Italia_merge$Label %in% TEs_in_gene_and_5kB_surroundings_A_Italia[,4],c(1:4,15,20)]
Retrotransposons_All_age_A_East_1kB_and_5kB <- All_A_East_merge[All_A_East_merge$Label %in% Retrotransposon_ID & !All_A_East_merge$Position %in% not_derived_TE_All_A_East[,2] & !All_A_East_merge$Label %in% TEs_in_gene_and_1kB_surroundings_A_East[,4] & All_A_East_merge$Label %in% TEs_in_gene_and_5kB_surroundings_A_East[,4],c(1:4,15,20)]
Retrotransposons_All_age_B_West_1kB_and_5kB <- All_B_West_merge[All_B_West_merge$Label %in% Retrotransposon_ID & !All_B_West_merge$Position %in% not_derived_TE_All_B_West[,2] & !All_B_West_merge$Label %in% TEs_in_gene_and_1kB_surroundings_B_West[,4] & All_B_West_merge$Label %in% TEs_in_gene_and_5kB_surroundings_B_West[,4],c(1:4,15,20)]

Retrotransposons_All_age_B_East_more_than_5kb <- All_B_East_merge[All_B_East_merge$Label %in% Retrotransposon_ID & !All_B_East_merge$Position %in% not_derived_TE_All_B_East[,2] & !All_B_East_merge$Label %in% TEs_in_gene_and_5kB_surroundings_B_East[,4],c(1:4,15,20)]
Retrotransposons_All_age_A_Italia_more_than_5kb <- All_A_Italia_merge[All_A_Italia_merge$Label %in% Retrotransposon_ID & !All_A_Italia_merge$Position %in% not_derived_TE_All_A_Italia[,2] & !All_A_Italia_merge$Label %in% TEs_in_gene_and_5kB_surroundings_A_Italia[,4],c(1:4,15,20)]
Retrotransposons_All_age_A_East_more_than_5kb <- All_A_East_merge[All_A_East_merge$Label %in% Retrotransposon_ID & !All_A_East_merge$Position %in% not_derived_TE_All_A_East[,2] & !All_A_East_merge$Label %in% TEs_in_gene_and_5kB_surroundings_A_East[,4],c(1:4,15,20)]
Retrotransposons_All_age_B_West_more_than_5kb <- All_B_West_merge[All_B_West_merge$Label %in% Retrotransposon_ID & !All_B_West_merge$Position %in% not_derived_TE_All_B_West[,2] & !All_B_West_merge$Label %in% TEs_in_gene_and_5kB_surroundings_B_West[,4],c(1:4,15,20)]


# downsapled neutral data 
resampled_synonymous_All_age_freq_matched_B_East_in_gene_and_1kB <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_East, Retrotransposons_All_age_B_East_in_gene_and_1kB)
resampled_synonymous_All_age_freq_matched_A_Italia_in_gene_and_1kB <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_Italia, Retrotransposons_All_age_A_Italia_in_gene_and_1kB)
resampled_synonymous_All_age_freq_matched_A_East_in_gene_and_1kB <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_East, Retrotransposons_All_age_A_East_in_gene_and_1kB)
resampled_synonymous_All_age_freq_matched_B_West_in_gene_and_1kB <- downsampling.neutral.SNP.function(synonymous_All_B_West_age, Retrotransposons_All_age_B_West_in_gene_and_1kB)

resampled_synonymous_All_age_freq_matched_B_East_1kB_and_5kB <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_East, Retrotransposons_All_age_B_East_1kB_and_5kB)
resampled_synonymous_All_age_freq_matched_A_Italia_1kB_and_5kB <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_Italia, Retrotransposons_All_age_A_Italia_1kB_and_5kB)
resampled_synonymous_All_age_freq_matched_A_East_1kB_and_5kB <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_East, Retrotransposons_All_age_A_East_1kB_and_5kB)
resampled_synonymous_All_age_freq_matched_B_West_1kB_and_5kB <- downsampling.neutral.SNP.function(synonymous_All_B_West_age, Retrotransposons_All_age_B_West_1kB_and_5kB)

resampled_synonymous_All_age_freq_matched_B_East_more_than_5kb <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_East, Retrotransposons_All_age_B_East_more_than_5kb)
resampled_synonymous_All_age_freq_matched_A_Italia_more_than_5kb <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_Italia, Retrotransposons_All_age_A_Italia_more_than_5kb)
resampled_synonymous_All_age_freq_matched_A_East_more_than_5kb <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_East, Retrotransposons_All_age_A_East_more_than_5kb)
resampled_synonymous_All_age_freq_matched_B_West_more_than_5kb <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_West, Retrotransposons_All_age_B_West_more_than_5kb)


## sort data 
Retrotransposons_All_age_sorted_B_East_in_gene_and_1kB <- Retrotransposons_All_age_B_East_in_gene_and_1kB[order(Retrotransposons_All_age_B_East_in_gene_and_1kB$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_B_East_in_gene_and_1kB <- resampled_synonymous_All_age_freq_matched_B_East_in_gene_and_1kB[order(resampled_synonymous_All_age_freq_matched_B_East_in_gene_and_1kB[,5]),]
Retrotransposons_All_age_sorted_A_Italia_in_gene_and_1kB <- Retrotransposons_All_age_A_Italia_in_gene_and_1kB[order(Retrotransposons_All_age_A_Italia_in_gene_and_1kB$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_A_Italia_in_gene_and_1kB <- resampled_synonymous_All_age_freq_matched_A_Italia_in_gene_and_1kB[order(resampled_synonymous_All_age_freq_matched_A_Italia_in_gene_and_1kB[,5]),]
Retrotransposons_All_age_sorted_A_East_in_gene_and_1kB <- Retrotransposons_All_age_A_East_in_gene_and_1kB[order(Retrotransposons_All_age_A_East_in_gene_and_1kB$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_A_East_in_gene_and_1kB <- resampled_synonymous_All_age_freq_matched_A_East_in_gene_and_1kB[order(resampled_synonymous_All_age_freq_matched_A_East_in_gene_and_1kB[,5]),]
Retrotransposons_All_age_sorted_B_West_in_gene_and_1kB <- Retrotransposons_All_age_B_West_in_gene_and_1kB[order(Retrotransposons_All_age_B_West_in_gene_and_1kB$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_B_West_in_gene_and_1kB <- resampled_synonymous_All_age_freq_matched_B_West_in_gene_and_1kB[order(resampled_synonymous_All_age_freq_matched_B_West_in_gene_and_1kB[,5]),]

Retrotransposons_All_age_sorted_B_East_1kB_and_5kB <- Retrotransposons_All_age_B_East_1kB_and_5kB[order(Retrotransposons_All_age_B_East_1kB_and_5kB$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_B_East_1kB_and_5kB <- resampled_synonymous_All_age_freq_matched_B_East_1kB_and_5kB[order(resampled_synonymous_All_age_freq_matched_B_East_1kB_and_5kB[,5]),]
Retrotransposons_All_age_sorted_A_Italia_1kB_and_5kB <- Retrotransposons_All_age_A_Italia_1kB_and_5kB[order(Retrotransposons_All_age_A_Italia_1kB_and_5kB$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_A_Italia_1kB_and_5kB <- resampled_synonymous_All_age_freq_matched_A_Italia_1kB_and_5kB[order(resampled_synonymous_All_age_freq_matched_A_Italia_1kB_and_5kB[,5]),]
Retrotransposons_All_age_sorted_A_East_1kB_and_5kB <- Retrotransposons_All_age_A_East_1kB_and_5kB[order(Retrotransposons_All_age_A_East_1kB_and_5kB$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_A_East_1kB_and_5kB <- resampled_synonymous_All_age_freq_matched_A_East_1kB_and_5kB[order(resampled_synonymous_All_age_freq_matched_A_East_1kB_and_5kB[,5]),]
Retrotransposons_All_age_sorted_B_West_1kB_and_5kB <- Retrotransposons_All_age_B_West_1kB_and_5kB[order(Retrotransposons_All_age_B_West_1kB_and_5kB$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_B_West_1kB_and_5kB <- resampled_synonymous_All_age_freq_matched_B_West_1kB_and_5kB[order(resampled_synonymous_All_age_freq_matched_B_West_1kB_and_5kB[,5]),]

Retrotransposons_All_age_sorted_B_East_more_than_5kb <- Retrotransposons_All_age_B_East_more_than_5kb[order(Retrotransposons_All_age_B_East_more_than_5kb$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_B_East_more_than_5kb <- resampled_synonymous_All_age_freq_matched_B_East_more_than_5kb[order(resampled_synonymous_All_age_freq_matched_B_East_more_than_5kb[,5]),]
Retrotransposons_All_age_sorted_A_Italia_more_than_5kb <- Retrotransposons_All_age_A_Italia_more_than_5kb[order(Retrotransposons_All_age_A_Italia_more_than_5kb$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_A_Italia_more_than_5kb <- resampled_synonymous_All_age_freq_matched_A_Italia_more_than_5kb[order(resampled_synonymous_All_age_freq_matched_A_Italia_more_than_5kb[,5]),]
Retrotransposons_All_age_sorted_A_East_more_than_5kb <- Retrotransposons_All_age_A_East_more_than_5kb[order(Retrotransposons_All_age_A_East_more_than_5kb$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_A_East_more_than_5kb <- resampled_synonymous_All_age_freq_matched_A_East_more_than_5kb[order(resampled_synonymous_All_age_freq_matched_A_East_more_than_5kb[,5]),]
Retrotransposons_All_age_sorted_B_West_more_than_5kb <- Retrotransposons_All_age_B_West_more_than_5kb[order(Retrotransposons_All_age_B_West_more_than_5kb$PostMode), ]
resampled_synonymous_All_age_sorted_freq_matched_B_West_more_than_5kb <- resampled_synonymous_All_age_freq_matched_B_West_more_than_5kb[order(resampled_synonymous_All_age_freq_matched_B_West_more_than_5kb[,5]),]


## get bins
All_ages_B_East_in_gene_and_1kB <- sort(c(Retrotransposons_All_age_sorted_B_East_in_gene_and_1kB$PostMode, resampled_synonymous_All_age_sorted_freq_matched_B_East_in_gene_and_1kB[,5]))
my_age_n_B_East_in_gene_and_1kB <- length(All_ages_B_East_in_gene_and_1kB)
my_bonderies_B_East_in_gene_and_1kB <- c(0)
for (n in 1:Nbin) {
  my_bonderies_B_East_in_gene_and_1kB <- c(my_bonderies_B_East_in_gene_and_1kB, All_ages_B_East_in_gene_and_1kB[floor(n*my_age_n_B_East_in_gene_and_1kB/Nbin)])
  
}
All_ages_A_Italia_in_gene_and_1kB <- sort(c(Retrotransposons_All_age_sorted_A_Italia_in_gene_and_1kB$PostMode, resampled_synonymous_All_age_sorted_freq_matched_A_Italia_in_gene_and_1kB[,5]))
my_age_n_A_Italia_in_gene_and_1kB <- length(All_ages_A_Italia_in_gene_and_1kB)
my_bonderies_A_Italia_in_gene_and_1kB <- c(0)
for (n in 1:Nbin) {
  my_bonderies_A_Italia_in_gene_and_1kB <- c(my_bonderies_A_Italia_in_gene_and_1kB, All_ages_A_Italia_in_gene_and_1kB[floor(n*my_age_n_A_Italia_in_gene_and_1kB/Nbin)])
  
}
All_ages_A_East_in_gene_and_1kB <- sort(c(Retrotransposons_All_age_sorted_A_East_in_gene_and_1kB$PostMode, resampled_synonymous_All_age_sorted_freq_matched_A_East_in_gene_and_1kB[,5]))
my_age_n_A_East_in_gene_and_1kB <- length(All_ages_A_East_in_gene_and_1kB)
my_bonderies_A_East_in_gene_and_1kB <- c(0)
for (n in 1:Nbin) {
  my_bonderies_A_East_in_gene_and_1kB <- c(my_bonderies_A_East_in_gene_and_1kB, All_ages_A_East_in_gene_and_1kB[floor(n*my_age_n_A_East_in_gene_and_1kB/Nbin)])
  
}
All_ages_B_West_in_gene_and_1kB <- sort(c(Retrotransposons_All_age_sorted_B_West_in_gene_and_1kB$PostMode, resampled_synonymous_All_age_sorted_freq_matched_B_West_in_gene_and_1kB[,5]))
my_age_n_B_West_in_gene_and_1kB <- length(All_ages_B_West_in_gene_and_1kB)
my_bonderies_B_West_in_gene_and_1kB <- c(0)
for (n in 1:Nbin) {
  my_bonderies_B_West_in_gene_and_1kB <- c(my_bonderies_B_West_in_gene_and_1kB, All_ages_B_West_in_gene_and_1kB[floor(n*my_age_n_B_West_in_gene_and_1kB/Nbin)])
  
}

All_ages_B_East_1kB_and_5kB <- sort(c(Retrotransposons_All_age_sorted_B_East_1kB_and_5kB$PostMode, resampled_synonymous_All_age_sorted_freq_matched_B_East_1kB_and_5kB[,5]))
my_age_n_B_East_1kB_and_5kB <- length(All_ages_B_East_1kB_and_5kB)
my_bonderies_B_East_1kB_and_5kB <- c(0)
for (n in 1:Nbin) {
  my_bonderies_B_East_1kB_and_5kB <- c(my_bonderies_B_East_1kB_and_5kB, All_ages_B_East_1kB_and_5kB[floor(n*my_age_n_B_East_1kB_and_5kB/Nbin)])
  
}
All_ages_A_Italia_1kB_and_5kB <- sort(c(Retrotransposons_All_age_sorted_A_Italia_1kB_and_5kB$PostMode, resampled_synonymous_All_age_sorted_freq_matched_A_Italia_1kB_and_5kB[,5]))
my_age_n_A_Italia_1kB_and_5kB <- length(All_ages_A_Italia_1kB_and_5kB)
my_bonderies_A_Italia_1kB_and_5kB <- c(0)
for (n in 1:Nbin) {
  my_bonderies_A_Italia_1kB_and_5kB <- c(my_bonderies_A_Italia_1kB_and_5kB, All_ages_A_Italia_1kB_and_5kB[floor(n*my_age_n_A_Italia_1kB_and_5kB/Nbin)])
  
}
All_ages_A_East_1kB_and_5kB <- sort(c(Retrotransposons_All_age_sorted_A_East_1kB_and_5kB$PostMode, resampled_synonymous_All_age_sorted_freq_matched_A_East_1kB_and_5kB[,5]))
my_age_n_A_East_1kB_and_5kB <- length(All_ages_A_East_1kB_and_5kB)
my_bonderies_A_East_1kB_and_5kB <- c(0)
for (n in 1:Nbin) {
  my_bonderies_A_East_1kB_and_5kB <- c(my_bonderies_A_East_1kB_and_5kB, All_ages_A_East_1kB_and_5kB[floor(n*my_age_n_A_East_1kB_and_5kB/Nbin)])
  
}
All_ages_B_West_1kB_and_5kB <- sort(c(Retrotransposons_All_age_sorted_B_West_1kB_and_5kB$PostMode, resampled_synonymous_All_age_sorted_freq_matched_B_West_1kB_and_5kB[,5]))
my_age_n_B_West_1kB_and_5kB <- length(All_ages_B_West_1kB_and_5kB)
my_bonderies_B_West_1kB_and_5kB <- c(0)
for (n in 1:Nbin) {
  my_bonderies_B_West_1kB_and_5kB <- c(my_bonderies_B_West_1kB_and_5kB, All_ages_B_West_1kB_and_5kB[floor(n*my_age_n_B_West_1kB_and_5kB/Nbin)])
  
}

All_ages_B_East_more_than_5kb <- sort(c(Retrotransposons_All_age_sorted_B_East_more_than_5kb$PostMode, resampled_synonymous_All_age_sorted_freq_matched_B_East_more_than_5kb[,5]))
my_age_n_B_East_more_than_5kb <- length(All_ages_B_East_more_than_5kb)
my_bonderies_B_East_more_than_5kb <- c(0)
for (n in 1:Nbin) {
  my_bonderies_B_East_more_than_5kb <- c(my_bonderies_B_East_more_than_5kb, All_ages_B_East_more_than_5kb[floor(n*my_age_n_B_East_more_than_5kb/Nbin)])
  
}
All_ages_A_Italia_more_than_5kb <- sort(c(Retrotransposons_All_age_sorted_A_Italia_more_than_5kb$PostMode, resampled_synonymous_All_age_sorted_freq_matched_A_Italia_more_than_5kb[,5]))
my_age_n_A_Italia_more_than_5kb <- length(All_ages_A_Italia_more_than_5kb)
my_bonderies_A_Italia_more_than_5kb <- c(0)
for (n in 1:Nbin) {
  my_bonderies_A_Italia_more_than_5kb <- c(my_bonderies_A_Italia_more_than_5kb, All_ages_A_Italia_more_than_5kb[floor(n*my_age_n_A_Italia_more_than_5kb/Nbin)])
  
}
All_ages_A_East_more_than_5kb <- sort(c(Retrotransposons_All_age_sorted_A_East_more_than_5kb$PostMode, resampled_synonymous_All_age_sorted_freq_matched_A_East_more_than_5kb[,5]))
my_age_n_A_East_more_than_5kb <- length(All_ages_A_East_more_than_5kb)
my_bonderies_A_East_more_than_5kb <- c(0)
for (n in 1:Nbin) {
  my_bonderies_A_East_more_than_5kb <- c(my_bonderies_A_East_more_than_5kb, All_ages_A_East_more_than_5kb[floor(n*my_age_n_A_East_more_than_5kb/Nbin)])
  
}
All_ages_B_West_more_than_5kb <- sort(c(Retrotransposons_All_age_sorted_B_West_more_than_5kb$PostMode, resampled_synonymous_All_age_sorted_freq_matched_B_West_more_than_5kb[,5]))
my_age_n_B_West_more_than_5kb <- length(All_ages_B_West_more_than_5kb)
my_bonderies_B_West_more_than_5kb <- c(0)
for (n in 1:Nbin) {
  my_bonderies_B_West_more_than_5kb <- c(my_bonderies_B_West_more_than_5kb, All_ages_B_West_more_than_5kb[floor(n*my_age_n_B_West_more_than_5kb/Nbin)])
  
}


## get delta frequency:
delta_freq_decile_All_resample_freq_matched_B_East_in_gene_and_1kB <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_B_East_in_gene_and_1kB, Retrotransposons_All_age_sorted_B_East_in_gene_and_1kB, my_bonderies_B_East_in_gene_and_1kB)
delta_freq_decile_All_resample_freq_matched_A_Italia_in_gene_and_1kB <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_A_Italia_in_gene_and_1kB, Retrotransposons_All_age_sorted_A_Italia_in_gene_and_1kB, my_bonderies_A_Italia_in_gene_and_1kB)
delta_freq_decile_All_resample_freq_matched_A_East_in_gene_and_1kB <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_A_East_in_gene_and_1kB, Retrotransposons_All_age_sorted_A_East_in_gene_and_1kB, my_bonderies_A_East_in_gene_and_1kB)
delta_freq_decile_All_resample_freq_matched_B_West_in_gene_and_1kB <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_B_West_in_gene_and_1kB, Retrotransposons_All_age_sorted_B_West_in_gene_and_1kB, my_bonderies_B_West_in_gene_and_1kB)

delta_freq_decile_All_resample_freq_matched_B_East_1kB_and_5kB <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_B_East_1kB_and_5kB, Retrotransposons_All_age_sorted_B_East_1kB_and_5kB, my_bonderies_B_East_1kB_and_5kB)
delta_freq_decile_All_resample_freq_matched_A_Italia_1kB_and_5kB <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_A_Italia_1kB_and_5kB, Retrotransposons_All_age_sorted_A_Italia_1kB_and_5kB, my_bonderies_A_Italia_1kB_and_5kB)
delta_freq_decile_All_resample_freq_matched_A_East_1kB_and_5kB <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_A_East_1kB_and_5kB, Retrotransposons_All_age_sorted_A_East_1kB_and_5kB, my_bonderies_A_East_1kB_and_5kB)
delta_freq_decile_All_resample_freq_matched_B_West_1kB_and_5kB <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_B_West_1kB_and_5kB, Retrotransposons_All_age_sorted_B_West_1kB_and_5kB, my_bonderies_B_West_1kB_and_5kB)

delta_freq_decile_All_resample_freq_matched_B_East_more_than_5kb <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_B_East_more_than_5kb, Retrotransposons_All_age_sorted_B_East_more_than_5kb, my_bonderies_B_East_more_than_5kb)
delta_freq_decile_All_resample_freq_matched_A_Italia_more_than_5kb <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_A_Italia_more_than_5kb, Retrotransposons_All_age_sorted_A_Italia_more_than_5kb, my_bonderies_A_Italia_more_than_5kb)
delta_freq_decile_All_resample_freq_matched_A_East_more_than_5kb <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_A_East_more_than_5kb, Retrotransposons_All_age_sorted_A_East_more_than_5kb, my_bonderies_A_East_more_than_5kb)
delta_freq_decile_All_resample_freq_matched_B_West_more_than_5kb <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_sorted_freq_matched_B_West_more_than_5kb, Retrotransposons_All_age_sorted_B_West_more_than_5kb, my_bonderies_B_West_more_than_5kb)


## bootstrapping
delta_freq_decile_All_boot_freq_matched_B_East_in_gene_and_1kB <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_B_East_in_gene_and_1kB[sample(1:dim(Retrotransposons_All_age_B_East_in_gene_and_1kB)[1], dim(Retrotransposons_All_age_B_East_in_gene_and_1kB)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_East, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_East_in_gene_and_1kB)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_A_Italia_in_gene_and_1kB <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_A_Italia_in_gene_and_1kB[sample(1:dim(Retrotransposons_All_age_A_Italia_in_gene_and_1kB)[1], dim(Retrotransposons_All_age_A_Italia_in_gene_and_1kB)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_Italia, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_Italia_in_gene_and_1kB)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_A_East_in_gene_and_1kB <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_A_East_in_gene_and_1kB[sample(1:dim(Retrotransposons_All_age_A_East_in_gene_and_1kB)[1], dim(Retrotransposons_All_age_A_East_in_gene_and_1kB)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_East, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_East_in_gene_and_1kB)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_B_West_in_gene_and_1kB <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_B_West_in_gene_and_1kB[sample(1:dim(Retrotransposons_All_age_B_West_in_gene_and_1kB)[1], dim(Retrotransposons_All_age_B_West_in_gene_and_1kB)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_West, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_West_in_gene_and_1kB)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)

delta_freq_decile_All_boot_freq_matched_B_East_1kB_and_5kB <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_B_East_1kB_and_5kB[sample(1:dim(Retrotransposons_All_age_B_East_1kB_and_5kB)[1], dim(Retrotransposons_All_age_B_East_1kB_and_5kB)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_East, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_East_1kB_and_5kB)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_A_Italia_1kB_and_5kB <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_A_Italia_1kB_and_5kB[sample(1:dim(Retrotransposons_All_age_A_Italia_1kB_and_5kB)[1], dim(Retrotransposons_All_age_A_Italia_1kB_and_5kB)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_Italia, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_Italia_1kB_and_5kB)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_A_East_1kB_and_5kB <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_A_East_1kB_and_5kB[sample(1:dim(Retrotransposons_All_age_A_East_1kB_and_5kB)[1], dim(Retrotransposons_All_age_A_East_1kB_and_5kB)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_East, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_East_1kB_and_5kB)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_B_West_1kB_and_5kB <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_B_West_1kB_and_5kB[sample(1:dim(Retrotransposons_All_age_B_West_1kB_and_5kB)[1], dim(Retrotransposons_All_age_B_West_1kB_and_5kB)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_West, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_West_1kB_and_5kB)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)

delta_freq_decile_All_boot_freq_matched_B_East_more_than_5kb <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_B_East_more_than_5kb[sample(1:dim(Retrotransposons_All_age_B_East_more_than_5kb)[1], dim(Retrotransposons_All_age_B_East_more_than_5kb)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_East, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_East_more_than_5kb)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_A_Italia_more_than_5kb <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_A_Italia_more_than_5kb[sample(1:dim(Retrotransposons_All_age_A_Italia_more_than_5kb)[1], dim(Retrotransposons_All_age_A_Italia_more_than_5kb)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_Italia, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_Italia_more_than_5kb)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_A_East_more_than_5kb <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_A_East_more_than_5kb[sample(1:dim(Retrotransposons_All_age_A_East_more_than_5kb)[1], dim(Retrotransposons_All_age_A_East_more_than_5kb)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_A_East, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_A_East_more_than_5kb)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)
delta_freq_decile_All_boot_freq_matched_B_West_more_than_5kb <- matrix(unlist(sapply(1:n_boot, function(x){
  Retrotransposons_All_age_resample <- Retrotransposons_All_age_B_West_more_than_5kb[sample(1:dim(Retrotransposons_All_age_B_West_more_than_5kb)[1], dim(Retrotransposons_All_age_B_West_more_than_5kb)[1], replace = TRUE),]
  resampled_synonymous_All_age_resample <- downsampling.neutral.SNP.function(freq_matched_synonymous_All_age_B_West, Retrotransposons_All_age_resample)
  Retrotransposons_All_age_resample_sorted <- Retrotransposons_All_age_resample[order(Retrotransposons_All_age_resample$PostMode), ]
  resampled_synonymous_All_age_resample_sorted <- resampled_synonymous_All_age_resample[order(resampled_synonymous_All_age_resample[,5]),]
  delta_freq_decile_All_resample <- decile.delta.age.fix.decile.size.variable.n.bin.function(resampled_synonymous_All_age_resample_sorted, Retrotransposons_All_age_resample_sorted, my_bonderies_B_West_more_than_5kb)
  print(x/n_boot)
  return(delta_freq_decile_All_resample)
})), nrow = n_boot, ncol = Nbin, byrow = TRUE)





