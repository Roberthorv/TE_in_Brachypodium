## get high iHS regions:

library(rehh)

## Read in iHS results from the result file with the computed iHS results:
B_West_iHS_Bd1 <- read.table("./scanB_WestIHS_Bd1.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_West_iHS_Bd1) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_West_iHS_Bd2 <- read.table("./scanB_WestIHS_Bd2.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_West_iHS_Bd2) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_West_iHS_Bd3 <- read.table("./scanB_WestIHS_Bd3.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_West_iHS_Bd3) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_West_iHS_Bd4 <- read.table("./scanB_WestIHS_Bd4.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_West_iHS_Bd4) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_West_iHS_Bd5 <- read.table("./scanB_WestIHS_Bd5.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_West_iHS_Bd5) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_West_iHS <- rbind(B_West_iHS_Bd1, B_West_iHS_Bd2, B_West_iHS_Bd3, B_West_iHS_Bd4, B_West_iHS_Bd5)

B_East_iHS_Bd1 <- read.table("./scanB_EastIHS_Bd1.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_East_iHS_Bd1) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_East_iHS_Bd2 <- read.table("./scanB_EastIHS_Bd2.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_East_iHS_Bd2) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_East_iHS_Bd3 <- read.table("./scanB_EastIHS_Bd3.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_East_iHS_Bd3) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_East_iHS_Bd4 <- read.table("./scanB_EastIHS_Bd4.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_East_iHS_Bd4) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_East_iHS_Bd5 <- read.table("./scanB_EastIHS_Bd5.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(B_East_iHS_Bd5) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
B_East_iHS <- rbind(B_East_iHS_Bd1, B_East_iHS_Bd2, B_East_iHS_Bd3, B_East_iHS_Bd4, B_East_iHS_Bd5)

A_East_iHS_Bd1 <- read.table("./scanA_EastIHS_Bd1.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_East_iHS_Bd1) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_East_iHS_Bd2 <- read.table("./scanA_EastIHS_Bd2.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_East_iHS_Bd2) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_East_iHS_Bd3 <- read.table("./scanA_EastIHS_Bd3.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_East_iHS_Bd3) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_East_iHS_Bd4 <- read.table("./scanA_EastIHS_Bd4.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_East_iHS_Bd4) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_East_iHS_Bd5 <- read.table("./scanA_EastIHS_Bd5.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_East_iHS_Bd5) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_East_iHS <- rbind(A_East_iHS_Bd1, A_East_iHS_Bd2, A_East_iHS_Bd3, A_East_iHS_Bd4, A_East_iHS_Bd5)

A_Italy_iHS_Bd1 <- read.table("./scanA_ItalyIHS_Bd1.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_Italy_iHS_Bd1) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_Italy_iHS_Bd2 <- read.table("./scanA_ItalyIHS_Bd2.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_Italy_iHS_Bd2) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_Italy_iHS_Bd3 <- read.table("./scanA_ItalyIHS_Bd3.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_Italy_iHS_Bd3) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_Italy_iHS_Bd4 <- read.table("./scanA_ItalyIHS_Bd4.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_Italy_iHS_Bd4) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_Italy_iHS_Bd5 <- read.table("./scanA_ItalyIHS_Bd5.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(A_Italy_iHS_Bd5) <- c("CHR", "POSITION", "IHS", "LOGPVALUE")
A_Italy_iHS <- rbind(A_Italy_iHS_Bd1, A_Italy_iHS_Bd2, A_Italy_iHS_Bd3, A_Italy_iHS_Bd4, A_Italy_iHS_Bd5)

## select the regions so that less than 5% of the genome is selected in each clade.
cr.B_East <- calc_candidate_regions(B_East_iHS,
                                    threshold = 4,
                                    pval = TRUE,
                                    window_size = 250000,
                                    overlap = 50000,
                                    min_n_extr_mrk = 10, 
                                    min_n_mrk = 10)

cr.B_West <- calc_candidate_regions(B_West_iHS,
                                    threshold = 6,
                                    pval = TRUE,
                                    window_size = 250000,
                                    overlap = 50000,
                                    min_n_extr_mrk = 10, 
                                    min_n_mrk = 10)

cr.A_East <- calc_candidate_regions(A_East_iHS,
                                    threshold = 6,
                                    pval = TRUE,
                                    window_size = 250000,
                                    overlap = 50000,
                                    min_n_extr_mrk = 10, 
                                    min_n_mrk = 10)

cr.A_Italy <- calc_candidate_regions(A_Italy_iHS,
                                     threshold = 5,
                                     pval = TRUE,
                                     window_size = 250000,
                                     overlap = 50000,
                                     min_n_extr_mrk = 10, 
                                     min_n_mrk = 10)

## Save iHS candidate regions in bed:
write.table(format(cr.B_West[,1:3], scientific = FALSE), file = "./high_iHS_regions_B_West_pos_sel_less_5_perc.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(format(cr.A_Italy[,1:3], scientific = FALSE), file = "./high_iHS_regions_A_Italy_pos_sel_less_5_perc.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(format(cr.A_East[,1:3], scientific = FALSE), file = "./high_iHS_regions_A_East_pos_sel_less_5_perc.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(format(cr.B_East[,1:3], scientific = FALSE), file = "./high_iHS_regions_B_East_pos_sel_less_5_perc.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)



