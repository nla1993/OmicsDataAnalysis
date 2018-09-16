# Occupancy of H1 variants at specific histone modifications 
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/PTMs/")

# H3K4me3
H1_2_H3K4me3<- read.table ('H3K4me3_GSM733682/H3K4me3H1_2_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1X_H3K4me3<- read.table ('H3K4me3_GSM733682/H3K4me3H1X_InputMergedSubtracted_fromBam.bed', header=FALSE)

H1.0_HA_H3K4me3<- read.table ('H3K4me3_GSM733682/H3K4me30_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1.2_HA_H3K4me3<- read.table ('H3K4me3_GSM733682/H3K4me32_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)


final_table <- cbind (H1_2_H3K4me3$V10, H1X_H3K4me3$V10, H1.0_HA_H3K4me3$V10,
                      H1.2_HA_H3K4me3$V10)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "H3K4me3", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                       "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)
#############################################
# H3K27me3
H1_2_H3K27me3<- read.table ('H3K27me3_GSM733696/H3K27me3H1_2_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1X_H3K27me3<- read.table ('H3K27me3_GSM733696/H3K27me3H1X_InputMergedSubtracted_fromBam.bed', header=FALSE)

H1.0_HA_H3K27me3<- read.table ('H3K27me3_GSM733696/H3K27me30_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1.2_HA_H3K27me3<- read.table ('H3K27me3_GSM733696/H3K27me32_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)

final_table <- cbind (H1_2_H3K27me3$V10, H1X_H3K27me3$V10, H1.0_HA_H3K27me3$V10,
                      H1.2_HA_H3K27me3$V10)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "H3K27me3", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                        "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)
#############################################
# H3K9me3
H1_2_H3K9me3<- read.table ('H3K9me3/H3K9me3H1_2_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1X_H3K9me3<- read.table ('H3K9me3/H3K9me3H1X_InputMergedSubtracted_fromBam.bed', header=FALSE)

H1.0_HA_H3K9me3<- read.table ('H3K9me3/H3K9me30_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1.2_HA_H3K9me3<- read.table ('H3K9me3/H3K9me32_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)

final_table <- cbind (H1_2_H3K9me3$V10, H1X_H3K9me3$V10, H1.0_HA_H3K9me3$V10,
                      H1.2_HA_H3K9me3$V10)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "H3K9me3", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                       "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)
#############################################
# H3K27ac
H1_2_H3K27ac<- read.table ('H3K27ac_GSM733684//H3K27acH1_2_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1X_H3K27ac<- read.table ('H3K27ac_GSM733684//H3K27acH1X_InputMergedSubtracted_fromBam.bed', header=FALSE)

H1.0_HA_H3K27ac<- read.table ('H3K27ac_GSM733684//H3K27ac0_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1.2_HA_H3K27ac<- read.table ('H3K27ac_GSM733684//H3K27ac2_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)

final_table <- cbind (H1_2_H3K27ac$V10, H1X_H3K27ac$V10, H1.0_HA_H3K27ac$V10,
                      H1.2_HA_H3K27ac$V10)
colnames(final_table)<- c("H1.2-endp", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "H3K27ac", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                       "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)
#############################################
# H3K36me3
H1_2_H3K36me3<- read.table ('H3K36me3_GSM733711/H3K36me3H1_2_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1X_H3K36me3<- read.table ('H3K36me3_GSM733711/H3K36me3H1X_InputMergedSubtracted_fromBam.bed', header=FALSE)

H1.0_HA_H3K36me3<- read.table ('H3K36me3_GSM733711/H3K36me30_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1.2_HA_H3K36me3<- read.table ('H3K36me3_GSM733711/H3K36me32_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)

final_table <- cbind (H1_2_H3K36me3$V10, H1X_H3K36me3$V10, H1.0_HA_H3K36me3$V10,
                      H1.2_HA_H3K36me3$V10)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "H3K36me3", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                        "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)
#############################################
# H3K4me1
H1_2_H3K4me1<- read.table ('H3K4me1_GSM798322/H3K4me1H1_2_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1X_H3K4me1<- read.table ('H3K4me1_GSM798322/H3K4me1H1X_InputMergedSubtracted_fromBam.bed', header=FALSE)

H1.0_HA_H3K4me1<- read.table ('H3K4me1_GSM798322/H3K4me10_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1.2_HA_H3K4me1<- read.table ('H3K4me1_GSM798322/H3K4me12_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)

final_table <- cbind (H1_2_H3K4me1$V10, H1X_H3K4me1$V10, H1.0_HA_H3K4me1$V10,
                      H1.2_HA_H3K4me1$V10)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "H3K4me1", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                       "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)
#############################################
# H3K9ac
H1_2_H3K9ac<- read.table ('H3K9ac_GSM733756/H3K9acH1_2_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1X_H3K9ac<- read.table ('H3K9ac_GSM733756/H3K9acH1X_InputMergedSubtracted_fromBam.bed', header=FALSE)

H1.0_HA_H3K9ac<- read.table ('H3K9ac_GSM733756/H3K9ac0_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)
H1.2_HA_H3K9ac<- read.table ('H3K9ac_GSM733756/H3K9ac2_HA_InputMergedSubtracted_fromBam.bed', header=FALSE)

final_table <- cbind (H1_2_H3K9ac$V10, H1X_H3K9ac$V10, H1.0_HA_H3K9ac$V10,
                      H1.2_HA_H3K9ac$V10)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "H3K9ac", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                      "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)

