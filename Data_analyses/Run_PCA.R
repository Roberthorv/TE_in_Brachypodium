## read in TIP TAP data:
Full_clean_TIP_TAP_Sample_table_342_samples <- read.table("./Full_clean_TIP_TAP_Sample_table_342_samples.txt", sep = "\t", header = TRUE)

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

## Total number of polymorphic TE insertions 
dim(Full_clean_TIP_TAP_Sample_table_326_samples)[1]

## Total number of retrotransposons and DNA transposons 
TIP_TAP_in_vcf_info_bed <- read.table("./TIP_TAP_in_vcf_info.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

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

Sample_info <- read.table("./Sample_info.txt")
Sample_info[Sample_info[,1] == "Bd21_Stritt", 2] <- "B_East"
Sample_info[Sample_info[,1] == "Bd21_Skalska", 2] <- "B_East"


## run PCA on SNPs

## read in data:
four_fold_degenerate_SNP_Bdis326_vcf <- read.vcfR("./four_fold_degenerate_SNP_Bdis326.vcf.gz", verbose = FALSE )

## extract SNPs
four_fold_degenerate_SNP_Bdis326_extract_gt <- extract.gt(four_fold_degenerate_SNP_Bdis326_vcf, element = "GT", IDtoRowNames  = F, as.numeric = T, convertNA = T, return.alleles = F)
four_fold_degenerate_SNP_Bdis326_extract_gt_df <- data.frame(t(four_fold_degenerate_SNP_Bdis326_extract_gt[complete.cases(four_fold_degenerate_SNP_Bdis326_extract_gt),]))

df_pca_data_four_fold_degenerate_SNP_Bdis326 <- four_fold_degenerate_SNP_Bdis326_extract_gt_df[,apply(four_fold_degenerate_SNP_Bdis326_extract_gt_df, 2, var, na.rm=TRUE) != 0]

pca_four_fold_degenerate_SNP_Bdis326.pca <- prcomp(df_pca_data_four_fold_degenerate_SNP_Bdis326, center = T, scale. = T)

Sample_info_b <- Sample_info

info_table <- read.table("./bdis_333_info.csv", sep = ",", header = TRUE)
info_table[info_table$source=="roulin",]

for (w in 1:dim(info_table[info_table$source=="roulin",])[1]) {
  Sample_info_b[Sample_info_b[,1] == info_table[info_table$source=="roulin",][w,1],1] <- info_table[info_table$source=="roulin",][w,5]
}

for (v in 1:dim(info_table[info_table$source=="danka",])[1]) {
  Sample_info_b[Sample_info_b[,1] == info_table[info_table$source=="danka",][v,1],1] <- info_table[info_table$source=="danka",][v,5]
}

data_326_samples_four_fold_degenerate_SNP.clade <- c()
for (i in 1:length(row.names(df_pca_data_four_fold_degenerate_SNP_Bdis326))) {
  data_326_samples_four_fold_degenerate_SNP.clade <- c(data_326_samples_four_fold_degenerate_SNP.clade, Sample_info_b[Sample_info_b$Names == gsub(pattern = "Bd3-1", replacement = "Bd3-1_r", row.names(df_pca_data_four_fold_degenerate_SNP_Bdis326))[i],]$Clade)
}

pdf(file = "./SNP_PCA.pdf")
ggbiplot(pca_four_fold_degenerate_SNP_Bdis326.pca, ellipse=TRUE,  labels=rownames(df_pca_data_four_fold_degenerate_SNP_Bdis326), groups=data_326_samples_four_fold_degenerate_SNP.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
dev.off()


## run PCA on TEs
pca_data_326_samples <- t(Full_clean_TIP_TAP_Sample_table_326_samples[complete.cases(Full_clean_TIP_TAP_Sample_table_326_samples[,11:336]), 11:336])
colnames(pca_data_326_samples) <- Full_clean_TIP_TAP_Sample_table_326_samples[complete.cases(Full_clean_TIP_TAP_Sample_table_326_samples[,11:336]),4]

df_pca_data_326_samples <- pca_data_326_samples[,apply(pca_data_326_samples, 2, var, na.rm=TRUE) != 0]
pca_326_samples.pca <- prcomp(df_pca_data_326_samples, center = T, scale. = T)

data_326_samples.clade <- c()
for (i in 1:length(row.names(df_pca_data_326_samples))) {
  data_326_samples.clade <- c(data_326_samples.clade, Sample_info[Sample_info$Names == gsub(pattern = "Bd3-1", replacement = "Bd3-1_r", gsub(pattern = "Bd21-3", replacement = "Bd21-3_r", gsub(pattern = "\\.", replacement = "-", gsub(pattern = "X", replacement = "", row.names(df_pca_data_326_samples)))))[i],]$Clade)
}

pdf(file = "./TE_PCA.pdf")
ggbiplot(pca_326_samples.pca,ellipse=TRUE,  labels=rownames(df_pca_data_326_samples), groups=data_326_samples.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
dev.off()

## only using TIPs
pca_TIP_data_326_samples <- pca_data_326_samples[,grep("TIP", colnames(pca_data_326_samples))]
df_pca_TIP_data_326_samples <- pca_TIP_data_326_samples[,apply(pca_TIP_data_326_samples, 2, var, na.rm=TRUE) != 0]
TIP_326_samples.pca <- prcomp(df_pca_TIP_data_326_samples, center = T, scale. = T)
pdf(file = "./TIP_PCA.pdf")
ggbiplot(TIP_326_samples.pca,ellipse=TRUE,  labels=rownames(df_pca_TIP_data_326_samples), groups=data_326_samples.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
dev.off()

## only using TAPs
pca_TAP_data_326_samples <- pca_data_326_samples[,grep("TAP", colnames(pca_data_326_samples))]
df_pca_TAP_data_326_samples <- pca_TAP_data_326_samples[,apply(pca_TAP_data_326_samples, 2, var, na.rm=TRUE) != 0]
TAP_326_samples.pca <- prcomp(df_pca_TAP_data_326_samples, center = T, scale. = T)
pdf(file = "./TAP_PCA.pdf")
ggbiplot(TAP_326_samples.pca,ellipse=TRUE,  labels=rownames(df_pca_TAP_data_326_samples), groups=data_326_samples.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
dev.off()

## plot SNP and TE PCA
TE_pca_plot <- ggbiplot(pca_326_samples.pca, ellipse=TRUE, groups=data_326_samples.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
SNP_pca_plot <- ggbiplot(pca_four_fold_degenerate_SNP_Bdis326.pca, ellipse=TRUE , groups=data_326_samples_four_fold_degenerate_SNP.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))

pdf(file = "./SNP_and_TE_PCA.pdf", width = 10, height = 3.5)
ggarrange(TE_pca_plot, SNP_pca_plot, 
          labels = c("Genetic variation in polymorphic TEs", "Genetic variation in synonymous SNPs"), label.x = c(-0.25,-0.25),
          ncol = 2, nrow = 1)
dev.off()


## run PCA on retroposons and DNA-transposons:

## only using Retrotransposon 
Full_clean_Retrotransposon_Sample_table_326_samples <- Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples$TE_ID %in% Retrotransposon_ID,]

pca_data_326_samples_Retrotransposon <- t(Full_clean_Retrotransposon_Sample_table_326_samples[complete.cases(Full_clean_Retrotransposon_Sample_table_326_samples[,11:336]), 11:336])
colnames(pca_data_326_samples_Retrotransposon) <- Full_clean_Retrotransposon_Sample_table_326_samples[complete.cases(Full_clean_Retrotransposon_Sample_table_326_samples[,11:336]),4]
df_pca_data_326_samples_Retrotransposon <- pca_data_326_samples_Retrotransposon[,apply(pca_data_326_samples_Retrotransposon, 2, var, na.rm=TRUE) != 0]
pca_326_samples_Retrotransposon.pca <- prcomp(df_pca_data_326_samples_Retrotransposon, center = T, scale. = T)

data_326_samples_Retrotransposon.clade <- c()
for (i in 1:length(row.names(df_pca_data_326_samples_Retrotransposon))) {
  data_326_samples_Retrotransposon.clade <- c(data_326_samples_Retrotransposon.clade, Sample_info[Sample_info$Names == gsub(pattern = "Bd3-1", replacement = "Bd3-1_r", gsub(pattern = "Bd21-3", replacement = "Bd21-3_r", gsub(pattern = "\\.", replacement = "-", gsub(pattern = "X", replacement = "", row.names(df_pca_data_326_samples_Retrotransposon)))))[i],]$Clade)
}

pdf(file = "./Retrotransposon_PCA.pdf")
ggbiplot(pca_326_samples_Retrotransposon.pca, ellipse=TRUE,  groups=data_326_samples_Retrotransposon.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
dev.off()

## only using DNAtransposon
Full_clean_DNAtransposon_Sample_table_326_samples <- Full_clean_TIP_TAP_Sample_table_326_samples[Full_clean_TIP_TAP_Sample_table_326_samples$TE_ID %in% DNAtransposon_ID,]

pca_data_326_samples_DNAtransposon <- t(Full_clean_DNAtransposon_Sample_table_326_samples[complete.cases(Full_clean_DNAtransposon_Sample_table_326_samples[,11:336]), 11:336])
colnames(pca_data_326_samples_DNAtransposon) <- Full_clean_DNAtransposon_Sample_table_326_samples[complete.cases(Full_clean_DNAtransposon_Sample_table_326_samples[,11:336]),4]
df_pca_data_326_samples_DNAtransposon <- pca_data_326_samples_DNAtransposon[,apply(pca_data_326_samples_DNAtransposon, 2, var, na.rm=TRUE) != 0]
pca_326_samples_DNAtransposon.pca <- prcomp(df_pca_data_326_samples_DNAtransposon, center = T, scale. = T)

data_326_samples_DNAtransposon.clade <- c()
for (i in 1:length(row.names(df_pca_data_326_samples_DNAtransposon))) {
  data_326_samples_DNAtransposon.clade <- c(data_326_samples_DNAtransposon.clade, Sample_info[Sample_info$Names == gsub(pattern = "Bd3-1", replacement = "Bd3-1_r", gsub(pattern = "Bd21-3", replacement = "Bd21-3_r", gsub(pattern = "\\.", replacement = "-", gsub(pattern = "X", replacement = "", row.names(df_pca_data_326_samples_DNAtransposon)))))[i],]$Clade)
}

pdf(file = "./DNAtransposon_PCA.pdf")
ggbiplot(pca_326_samples_DNAtransposon.pca, ellipse=TRUE, groups=data_326_samples_DNAtransposon.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
dev.off()


## plot SNP, TE, Retrotransposon and DNAtransposon PCA
TE_pca_plot <- ggbiplot(pca_326_samples.pca, ellipse=TRUE, groups=data_326_samples.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
SNP_pca_plot <- ggbiplot(pca_four_fold_degenerate_SNP_Bdis326.pca, ellipse=TRUE , groups=data_326_samples_four_fold_degenerate_SNP.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
Retrotransposon_pca_plot <- ggbiplot(pca_326_samples_Retrotransposon.pca, ellipse=TRUE,  groups=data_326_samples_Retrotransposon.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))
DNAtransposon_pca_plot <- ggbiplot(pca_326_samples_DNAtransposon.pca, ellipse=TRUE, groups=data_326_samples_DNAtransposon.clade, var.axes = F) +
  scale_color_manual(name="Clades", values = c("#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF"))

pdf(file = "./SNP_TE_Retrotransposon_and_DNAtransposon_PCA.pdf", width = 14, height = 11)
ggarrange(TE_pca_plot, SNP_pca_plot, Retrotransposon_pca_plot, DNAtransposon_pca_plot, 
          labels = c("Genetic variation in polymorphic TEs", "Genetic variation in synonymous SNPs", "Genetic variation in Retrotransposons", "Genetic variation in DNA-transposons"),
          ncol = 2, nrow = 2)
dev.off()

## plot Retrotransposon and DNAtransposon PCA
pdf(file = "./Retrotransposon_and_DNAtransposon_PCA.pdf", width = 10, height = 6.5)
ggarrange(Retrotransposon_pca_plot, DNAtransposon_pca_plot, 
          labels = c("Genetic variation in Retrotransposons", "Genetic variation in DNA-transposons"), label.x = c(-0.26,-0.24),
          ncol = 2, nrow = 1)
dev.off()


## mentel test:
library(ade4)

## DNA transposon
Mentel_test_DNAtransposon_data <- cbind(data_326_samples_DNAtransposon.clade, pca_326_samples_DNAtransposon.pca$x[,1:2])

DNAtransposon_PC.dists <- dist(cbind(as.numeric(Mentel_test_DNAtransposon_data[,2]), as.numeric(Mentel_test_DNAtransposon_data[,3])))
DNAtransposon_clade.dists <- dist(as.numeric(as.factor(Mentel_test_DNAtransposon_data[,1])))

mantel.rtest(DNAtransposon_PC.dists, DNAtransposon_clade.dists, nrepet = 99)

## Retrotransposon
Mentel_test_Retrotransposon_data <- cbind(data_326_samples_Retrotransposon.clade, pca_326_samples_Retrotransposon.pca$x[,1:2])

Retrotransposon_PC.dists <- dist(cbind(as.numeric(Mentel_test_Retrotransposon_data[,2]), as.numeric(Mentel_test_Retrotransposon_data[,3])))
Retrotransposon_clade.dists <- dist(as.numeric(as.factor(Mentel_test_Retrotransposon_data[,1])))

mantel.rtest(Retrotransposon_PC.dists, Retrotransposon_clade.dists, nrepet = 999)

