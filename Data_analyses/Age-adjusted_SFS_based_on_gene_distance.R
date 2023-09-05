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

## plotting the results 

my_names_A_East <- rep(NA, 45)
for (a in 1:45) {
  my_names_A_East[a] <- paste0(format(my_bonderies_A_East[floor((a-1)/3)+1], scientific = TRUE, digits = 2), " -\n", format(my_bonderies_A_East[floor((a-1)/3)+2], scientific = TRUE, digits = 2))
}
my_names_A_Italia <- rep(NA, 45)
for (a in 1:45) {
  my_names_A_Italia[a] <- paste0(format(my_bonderies_A_Italia[floor((a-1)/3)+1], scientific = TRUE, digits = 2), " -\n", format(my_bonderies_A_Italia[floor((a-1)/3)+2], scientific = TRUE, digits = 2))
}
my_names_B_West <- rep(NA, 45)
for (a in 1:45) {
  my_names_B_West[a] <- paste0(format(my_bonderies_B_West[floor((a-1)/3)+1], scientific = TRUE, digits = 2), " -\n", format(my_bonderies_B_West[floor((a-1)/3)+2], scientific = TRUE, digits = 2))
}
my_names_B_East <- rep(NA, 45)
for (a in 1:45) {
  my_names_B_East[a] <- paste0(format(my_bonderies_B_East[floor((a-1)/3)+1], scientific = TRUE, digits = 2), " -\n", format(my_bonderies_B_East[floor((a-1)/3)+2], scientific = TRUE, digits = 2))
}

library(scales)
my_at <- as.vector(rbind(1:Nbin, 1:Nbin, 1:Nbin))
my_plot_seq <- as.vector(rbind(c(1:Nbin), c((Nbin+1):(2*Nbin)), c((2*Nbin+1):(3*Nbin))))

pdf(file = "Delta_freqeuncy_freq_matched_TEs_vs_SNPs_combined.pdf", width = 11, height = 5)
layout(matrix(c(rep(1,19),rep(2,18),rep(3,18),rep(4,18),
                rep(1,19),rep(2,18),rep(3,18),rep(4,18),
                rep(5,19),rep(6,18),rep(7,18),rep(8,18),
                rep(5,19),rep(6,18),rep(7,18),rep(8,18),
                rep(5,19),rep(6,18),rep(7,18),rep(8,18)), nrow = 5, ncol = 18*4+1, byrow = TRUE))
par(mar=c(1, 4, 1.75, 0))

boxplot(cbind(delta_freq_decile_All_boot_A_East_non_synonymous, delta_freq_decile_All_boot_A_East_high_effect_SNPs, delta_freq_decile_All_boot_freq_matched_A_East)[,my_plot_seq],
        col = c(alpha("azure3", 0.7), alpha("azure4", 0.7), alpha("#96127d", 0.7)), xlab = "", ylab = "", cex.main = 2, main = "", 
        border = c("azure3", "azure4", "#96127d"),
        names = rep(NA, 45), cex.axis = 1.4,
        at = my_at, outline = FALSE, ylim = c(-0.45, 0.15))
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("Retrotransposons ","Non-synonymous SNPs","High effect SNPs"), 
       col = c("#96127d", "azure3", "azure4"), bty = "n", pch = 15, cex = 1.1)
title(ylab = expression(paste(Delta, " frequency")), mgp = c(2.3, 1, 0), cex.lab = 1.8)
mtext("A", side = 3, at = -1.5, line = -0.5, cex = 1.3)
mtext("***", side = 3, at = 12, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 13, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 14, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 15, line = -1.4, cex = 0.6)

par(mar=c(1, 2.5, 1.75, 0.5))

boxplot(cbind(delta_freq_decile_All_boot_A_Italia_non_synonymous, delta_freq_decile_All_boot_A_Italia_high_effect_SNPs, delta_freq_decile_All_boot_freq_matched_A_Italia)[,my_plot_seq],
        col = c(alpha("azure3", 0.7), alpha("azure4", 0.7), alpha("#21908c", 0.7)), xlab = "", ylab = "", cex.main = 2, main = "", yaxt = "n", 
        border = c("azure3", "azure4", "#21908c"), 
        names = rep(NA, 45), 
        at = my_at, outline=FALSE, ylim = c(-0.45, 0.15))
axis(2, col.axis = "white")
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("Retrotransposons ","Non-synonymous SNPs","High effect SNPs"), 
       col = c("#21908c", "azure3", "azure4"), bty = "n", pch = 15, cex = 1.1)
mtext("B", side = 3, at = -1, line = -0.5, cex = 1.3)
mtext("***", side = 3, at = 13, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 14, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 15, line = -1.4, cex = 0.6)

boxplot(cbind(delta_freq_decile_All_boot_B_West_non_synonymous, delta_freq_decile_All_boot_B_West_high_effect_SNPs, delta_freq_decile_All_boot_freq_matched_B_West)[,my_plot_seq],
        col = c(alpha("azure3", 0.7), alpha("azure4", 0.7), alpha("#3B528BFF", 0.7)), xlab = "", ylab = "", cex.main = 2, main = "", yaxt = "n",
        border = c("azure3", "azure4", "#3B528BFF"),
        names = rep(NA, 45), 
        at = my_at, outline=FALSE, ylim = c(-0.45, 0.15))
axis(2, col.axis = "white")
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("Retrotransposons ","Non-synonymous SNPs","High effect SNPs"), 
       col = c("#3B528BFF", "azure3", "azure4"), bty = "n", pch = 15, cex = 1.1)
mtext("C", side = 3, at = -1.2, line = -0.5, cex = 1.3)
mtext("***", side = 3, at = 11, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 12, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 13, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 14, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 15, line = -1.4, cex = 0.6)

boxplot(cbind(delta_freq_decile_All_boot_B_East_non_synonymous, delta_freq_decile_All_boot_B_East_high_effect_SNPs, delta_freq_decile_All_boot_freq_matched_B_East)[,my_plot_seq],
        col = c(alpha("azure3", 0.7), alpha("azure4", 0.7), alpha("#5DC863FF", 0.7)), xlab = "", ylab = "", cex.main = 2, main = "", yaxt = "n",
        border = c("azure3", "azure4", "#5DC863FF"),
        names = rep(NA, 45), 
        at = my_at, outline=FALSE, ylim = c(-0.45, 0.15))
axis(2, col.axis = "white")
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("Retrotransposons ","Non-synonymous SNPs","High effect SNPs"), 
       col = c("#5DC863FF", "azure3", "azure4"), bty = "n", pch = 15, cex = 1.1)
mtext("D", side = 3, at = -1.2, line = -0.5, cex = 1.3)
mtext("***", side = 3, at = 11, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 12, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 13, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 14, line = -1.4, cex = 0.6)
mtext("***", side = 3, at = 15, line = -1.4, cex = 0.6)

par(mar=c(9.5, 4, 1, 0))

boxplot(cbind(delta_freq_decile_All_boot_freq_matched_A_East_in_gene_and_1kB, delta_freq_decile_All_boot_freq_matched_A_East_1kB_and_5kB, delta_freq_decile_All_boot_freq_matched_A_East_more_than_5kb)[,my_plot_seq],
        col = c(alpha("#96127d", 1), alpha("#96127d", 0.75), alpha("#96127d", 0.5)), xlab = "", ylab = "", cex.main = 2, main = "", 
        border = c("black", "azure4", "azure3"),
        names = rep(NA, 45), cex.axis = 1.4,
        at = my_at, outline=FALSE, ylim = c(-0.45, 0.15))
title(xlab = "Age bins", mgp = c(8, 1, 0), cex.lab = 1.75)
title(ylab = expression(paste(Delta, " frequency")), mgp = c(2.3, 1, 0), cex.lab = 1.8)
text(x = my_at + 0.55, y = rep(-0.52, 45), my_names_A_East, srt = 75, xpd = TRUE, cex = 1.25, pos = 2)
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("In genes and 1 kb surrounding","Between 1 and 5 kb away from genes","More than 5 kb away from genes"), 
       fill = c(alpha("#96127d", 1), alpha("#96127d", 0.75), alpha("#96127d", 0.5)), bty = "n", border = c("black", "azure4", "azure3"), cex = 1)
mtext("***", side = 3, at = 9, line = -1.1, cex = 0.6, col = alpha("#96127d", 1))
mtext("***", side = 3, at = 10, line = -1.1, cex = 0.6, col = alpha("#96127d", 1))
mtext("***", side = 3, at = 12, line = -1.1, cex = 0.6, col = alpha("#96127d", 1))
mtext("***", side = 3, at = 13, line = -1.1, cex = 0.6, col = alpha("#96127d", 1))
mtext("***", side = 3, at = 14, line = -1.1, cex = 0.6, col = alpha("#96127d", 1))
mtext("***", side = 3, at = 15, line = -1.1, cex = 0.6, col = alpha("#96127d", 1))

mtext("***", side = 3, at = 12, line = -1.4, cex = 0.6, col = alpha("#96127d", 0.75))
mtext("***", side = 3, at = 13, line = -1.4, cex = 0.6, col = alpha("#96127d", 0.75))
mtext("***", side = 3, at = 14, line = -1.4, cex = 0.6, col = alpha("#96127d", 0.75))
mtext("***", side = 3, at = 15, line = -1.4, cex = 0.6, col = alpha("#96127d", 0.75))

mtext("***", side = 3, at = 13, line = -1.7, cex = 0.6, col = alpha("#96127d", 0.5))
mtext("***", side = 3, at = 14, line = -1.7, cex = 0.6, col = alpha("#96127d", 0.5))
mtext("***", side = 3, at = 15, line = -1.7, cex = 0.6, col = alpha("#96127d", 0.5))

par(mar=c(9.5, 2.5, 1, 0.5))

boxplot(cbind(delta_freq_decile_All_boot_freq_matched_A_Italia_in_gene_and_1kB, delta_freq_decile_All_boot_freq_matched_A_Italia_1kB_and_5kB, delta_freq_decile_All_boot_freq_matched_A_Italia_more_than_5kb)[,my_plot_seq],
        col = c(alpha("#21908c", 1), alpha("#21908c", 0.75), alpha("#21908c", 0.5)), xlab = "", ylab = "", cex.main = 2, main = "",  yaxt = "n",
        border = c("black", "azure4", "azure3"), 
        names = rep(NA, 45), 
        at = my_at, outline=FALSE, ylim = c(-0.45, 0.15))
axis(2, col.axis = "white")
title(xlab = "Age bins", mgp = c(8, 1, 0), cex.lab = 1.75)
text(x = my_at + 0.55, y = rep(-0.52, 45), my_names_A_Italia, srt = 75, xpd = TRUE, cex = 1.25, pos = 2)
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("In genes and 1 kb surrounding","Between 1 and 5 kb away from genes","More than 5 kb away from genes"), 
       fill = c(alpha("#21908c", 1), alpha("#21908c", 0.75), alpha("#21908c", 0.5)), bty = "n", border = c("black", "azure4", "azure3"), cex = 1)
mtext("***", side = 3, at = 2, line = -1.1, cex = 0.6, col = alpha("#21908c", 1))
mtext("***", side = 3, at = 9, line = -1.1, cex = 0.6, col = alpha("#21908c", 1))
mtext("***", side = 3, at = 12, line = -1.1, cex = 0.6, col = alpha("#21908c", 1))
mtext("***", side = 3, at = 13, line = -1.1, cex = 0.6, col = alpha("#21908c", 1))
mtext("***", side = 3, at = 14, line = -1.1, cex = 0.6, col = alpha("#21908c", 1))
mtext("***", side = 3, at = 15, line = -1.1, cex = 0.6, col = alpha("#21908c", 1))

mtext("***", side = 3, at = 2, line = -1.4, cex = 0.6, col = alpha("#21908c", 0.75))
mtext("***", side = 3, at = 12, line = -1.4, cex = 0.6, col = alpha("#21908c", 0.75))
mtext("***", side = 3, at = 14, line = -1.4, cex = 0.6, col = alpha("#21908c", 0.75))
mtext("***", side = 3, at = 15, line = -1.4, cex = 0.6, col = alpha("#21908c", 0.75))

mtext("***", side = 3, at = 12, line = -1.7, cex = 0.6, col = alpha("#21908c", 0.5))
mtext("***", side = 3, at = 14, line = -1.7, cex = 0.6, col = alpha("#21908c", 0.5))
mtext("***", side = 3, at = 15, line = -1.7, cex = 0.6, col = alpha("#21908c", 0.5))


boxplot(cbind(delta_freq_decile_All_boot_freq_matched_B_West_in_gene_and_1kB, delta_freq_decile_All_boot_freq_matched_B_West_1kB_and_5kB, delta_freq_decile_All_boot_freq_matched_B_West_more_than_5kb)[,my_plot_seq],
        col = c(alpha("#3B528BFF", 1), alpha("#3B528BFF", 0.75), alpha("#3B528BFF", 0.5)), xlab = "", ylab = "", cex.main = 2, main = "",  yaxt = "n",
        border = c("black", "azure4", "azure3"),
        names = rep(NA, 45),
        at = my_at, outline=FALSE, ylim = c(-0.45, 0.15))
axis(2, col.axis = "white")
title(xlab = "Age bins", mgp = c(8, 1, 0), cex.lab = 1.75)
text(x = my_at + 0.55, y = rep(-0.52, 45), my_names_B_West, srt = 75, xpd = TRUE, cex = 1.25, pos = 2)
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("In genes and 1 kb surrounding","Between 1 and 5 kb away from genes","More than 5 kb away from genes"), 
       fill = c(alpha("#3B528BFF", 1), alpha("#3B528BFF", 0.75), alpha("#3B528BFF", 0.5)), bty = "n", border = c("black", "azure4", "azure3"), cex = 1)
mtext("***", side = 3, at = 10, line = -1.1, cex = 0.6, col = alpha("#3B528BFF", 1))
mtext("***", side = 3, at = 11, line = -1.1, cex = 0.6, col = alpha("#3B528BFF", 1))
mtext("***", side = 3, at = 12, line = -1.1, cex = 0.6, col = alpha("#3B528BFF", 1))
mtext("***", side = 3, at = 13, line = -1.1, cex = 0.6, col = alpha("#3B528BFF", 1))
mtext("***", side = 3, at = 14, line = -1.1, cex = 0.6, col = alpha("#3B528BFF", 1))
mtext("***", side = 3, at = 15, line = -1.1, cex = 0.6, col = alpha("#3B528BFF", 1))

mtext("***", side = 3, at = 2, line = -1.4, cex = 0.6, col = alpha("#3B528BFF", 0.75))
mtext("***", side = 3, at = 10, line = -1.4, cex = 0.6, col = alpha("#3B528BFF", 0.75))
mtext("***", side = 3, at = 11, line = -1.4, cex = 0.6, col = alpha("#3B528BFF", 0.75))
mtext("***", side = 3, at = 12, line = -1.4, cex = 0.6, col = alpha("#3B528BFF", 0.75))
mtext("***", side = 3, at = 13, line = -1.4, cex = 0.6, col = alpha("#3B528BFF", 0.75))
mtext("***", side = 3, at = 14, line = -1.4, cex = 0.6, col = alpha("#3B528BFF", 0.75))
mtext("***", side = 3, at = 15, line = -1.4, cex = 0.6, col = alpha("#3B528BFF", 0.75))

mtext("***", side = 3, at = 1, line = -1.7, cex = 0.6, col = alpha("#3B528BFF", 0.5))
mtext("***", side = 3, at = 3, line = -1.7, cex = 0.6, col = alpha("#3B528BFF", 0.5))
mtext("***", side = 3, at = 14, line = -1.7, cex = 0.6, col = alpha("#3B528BFF", 0.5))
mtext("***", side = 3, at = 15, line = -1.7, cex = 0.6, col = alpha("#3B528BFF", 0.5))


boxplot(cbind(delta_freq_decile_All_boot_freq_matched_B_East_in_gene_and_1kB, delta_freq_decile_All_boot_freq_matched_B_East_1kB_and_5kB, delta_freq_decile_All_boot_freq_matched_B_East_more_than_5kb)[,my_plot_seq],
        col = c(alpha("#5DC863FF", 1), alpha("#5DC863FF", 0.75), alpha("#5DC863FF", 0.5)), xlab = "", ylab = "", cex.main = 2, main = "",  yaxt = "n",
        border = c("black", "azure4", "azure3"),
        names = rep(NA, 45), 
        at = my_at, outline=FALSE, ylim = c(-0.45, 0.15))
axis(2, col.axis = "white")
title(xlab = "Age bins", mgp = c(8, 1, 0), cex.lab = 1.75)
text(x = my_at + 0.55, y = rep(-0.52, 45), my_names_B_East, srt = 75, xpd = TRUE, cex = 1.25, pos = 2)
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("In genes and 1 kb surrounding","Between 1 and 5 kb away from genes","More than 5 kb away from genes"), 
       fill = c(alpha("#5DC863FF", 1), alpha("#5DC863FF", 0.75), alpha("#5DC863FF", 0.5)), bty = "n", border = c("black", "azure4", "azure3"), cex = 1)
mtext("***", side = 3, at = 12, line = -1.1, cex = 0.6, col = alpha("#5DC863FF", 1))
mtext("***", side = 3, at = 13, line = -1.1, cex = 0.6, col = alpha("#5DC863FF", 1))
mtext("***", side = 3, at = 14, line = -1.1, cex = 0.6, col = alpha("#5DC863FF", 1))
mtext("***", side = 3, at = 15, line = -1.1, cex = 0.6, col = alpha("#5DC863FF", 1))

mtext("***", side = 3, at = 2, line = -1.4, cex = 0.6, col = alpha("#5DC863FF", 0.75))
mtext("***", side = 3, at = 9, line = -1.4, cex = 0.6, col = alpha("#5DC863FF", 0.75))
mtext("***", side = 3, at = 11, line = -1.4, cex = 0.6, col = alpha("#5DC863FF", 0.75))
mtext("***", side = 3, at = 12, line = -1.4, cex = 0.6, col = alpha("#5DC863FF", 0.75))
mtext("***", side = 3, at = 13, line = -1.4, cex = 0.6, col = alpha("#5DC863FF", 0.75))
mtext("***", side = 3, at = 14, line = -1.4, cex = 0.6, col = alpha("#5DC863FF", 0.75))
mtext("***", side = 3, at = 15, line = -1.4, cex = 0.6, col = alpha("#5DC863FF", 0.75))

mtext("***", side = 3, at = 1, line = -1.7, cex = 0.6, col = alpha("#5DC863FF", 0.5))
mtext("***", side = 3, at = 2, line = -1.7, cex = 0.6, col = alpha("#5DC863FF", 0.5))
mtext("***", side = 3, at = 9, line = -1.7, cex = 0.6, col = alpha("#5DC863FF", 0.5))
mtext("***", side = 3, at = 11, line = -1.7, cex = 0.6, col = alpha("#5DC863FF", 0.5))
mtext("***", side = 3, at = 13, line = -1.7, cex = 0.6, col = alpha("#5DC863FF", 0.5))
mtext("***", side = 3, at = 14, line = -1.7, cex = 0.6, col = alpha("#5DC863FF", 0.5))
mtext("***", side = 3, at = 15, line = -1.7, cex = 0.6, col = alpha("#5DC863FF", 0.5))

dev.off()





