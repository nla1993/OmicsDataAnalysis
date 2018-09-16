# PTMs and CRs data
# Normalized by length
# Reading the files


setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/normalized_by_cytoband_length/")
library(pheatmap)
library(RColorBrewer)

CHD1_gpos25_FDR <- read.table ('CHD1_FDR_0.01cytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
CHD1_gpos25_FDR$V8 <- paste(CHD1_gpos25_FDR$V4, CHD1_gpos25_FDR$V5, CHD1_gpos25_FDR$V1, sep="_")
CHD1_gpos50_FDR <- read.table ('CHD1_FDR_0.01cytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
CHD1_gpos50_FDR$V8 <- paste(CHD1_gpos50_FDR$V4, CHD1_gpos50_FDR$V5, CHD1_gpos50_FDR$V1, sep="_")
CHD1_gpos75_FDR <- read.table ('CHD1_FDR_0.01cytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
CHD1_gpos75_FDR$V8 <- paste(CHD1_gpos75_FDR$V4, CHD1_gpos75_FDR$V5, CHD1_gpos75_FDR$V1, sep="_")
CHD1_gpos100_FDR <- read.table ('CHD1_FDR_0.01cytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
CHD1_gpos100_FDR$V8 <- paste(CHD1_gpos100_FDR$V4, CHD1_gpos100_FDR$V5, CHD1_gpos100_FDR$V1, sep="_")


P300_gpos25 <- read.table ('P300_ENCFF001VIZcytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
P300_gpos25$V8 <- paste(P300_gpos25$V4, P300_gpos25$V5, P300_gpos25$V1, sep="_")
P300_gpos50 <- read.table ('P300_ENCFF001VIZcytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
P300_gpos50$V8 <- paste(P300_gpos50$V4, P300_gpos50$V5, P300_gpos50$V1, sep="_")
P300_gpos75 <- read.table ('P300_ENCFF001VIZcytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
P300_gpos75$V8 <- paste(P300_gpos75$V4, P300_gpos75$V5, P300_gpos75$V1, sep="_")
P300_gpos100 <- read.table ('P300_ENCFF001VIZcytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
P300_gpos100$V8 <- paste(P300_gpos100$V4, P300_gpos100$V5, P300_gpos100$V1, sep="_")


EZH2_gpos25 <- read.table ('EZH2_GSM1003520cytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
EZH2_gpos25$V8 <- paste(EZH2_gpos25$V4, EZH2_gpos25$V5, EZH2_gpos25$V1, sep="_")
EZH2_gpos50 <- read.table ('EZH2_GSM1003520cytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
EZH2_gpos50$V8 <- paste(EZH2_gpos50$V4, EZH2_gpos50$V5, EZH2_gpos50$V1, sep="_")
EZH2_gpos75 <- read.table ('EZH2_GSM1003520cytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
EZH2_gpos75$V8 <- paste(EZH2_gpos75$V4, EZH2_gpos75$V5, EZH2_gpos75$V1, sep="_")
EZH2_gpos100 <- read.table ('EZH2_GSM1003520cytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
EZH2_gpos100$V8 <- paste(EZH2_gpos100$V4, EZH2_gpos100$V5, EZH2_gpos100$V1, sep="_")

H3K27ac_gpos25 <- read.table ('H3K27accytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
H3K27ac_gpos25$V8 <- paste(H3K27ac_gpos25$V4, H3K27ac_gpos25$V5, H3K27ac_gpos25$V1, sep="_")
H3K27ac_gpos50 <- read.table ('H3K27accytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
H3K27ac_gpos50$V8 <- paste(H3K27ac_gpos50$V4, H3K27ac_gpos50$V5, H3K27ac_gpos50$V1, sep="_")
H3K27ac_gpos75 <- read.table ('H3K27accytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
H3K27ac_gpos75$V8 <- paste(H3K27ac_gpos75$V4, H3K27ac_gpos75$V5, H3K27ac_gpos75$V1, sep="_")
H3K27ac_gpos100 <- read.table ('H3K27accytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
H3K27ac_gpos100$V8 <- paste(H3K27ac_gpos100$V4, H3K27ac_gpos100$V5, H3K27ac_gpos100$V1, sep="_")

H3K27me3_gpos25 <- read.table ('H3K27me3_GSM733696cytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
H3K27me3_gpos25$V8 <- paste(H3K27me3_gpos25$V4, H3K27me3_gpos25$V5, H3K27me3_gpos25$V1, sep="_")
H3K27me3_gpos50 <- read.table ('H3K27me3_GSM733696cytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
H3K27me3_gpos50$V8 <- paste(H3K27me3_gpos50$V4, H3K27me3_gpos50$V5, H3K27me3_gpos50$V1, sep="_")
H3K27me3_gpos75 <- read.table ('H3K27me3_GSM733696cytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
H3K27me3_gpos75$V8 <- paste(H3K27me3_gpos75$V4, H3K27me3_gpos75$V5, H3K27me3_gpos75$V1, sep="_")
H3K27me3_gpos100 <- read.table ('H3K27me3_GSM733696cytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
H3K27me3_gpos100$V8 <- paste(H3K27me3_gpos100$V4, H3K27me3_gpos100$V5, H3K27me3_gpos100$V1, sep="_")

H3K4me3_gpos25 <- read.table ('H3K4me3_GSM733682cytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
H3K4me3_gpos25$V8 <- paste(H3K4me3_gpos25$V4, H3K4me3_gpos25$V5, H3K4me3_gpos25$V1, sep="_")
H3K4me3_gpos50 <- read.table ('H3K4me3_GSM733682cytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
H3K4me3_gpos50$V8 <- paste(H3K4me3_gpos50$V4, H3K4me3_gpos50$V5, H3K4me3_gpos50$V1, sep="_")
H3K4me3_gpos75 <- read.table ('H3K4me3_GSM733682cytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
H3K4me3_gpos75$V8 <- paste(H3K4me3_gpos75$V4, H3K4me3_gpos75$V5, H3K4me3_gpos75$V1, sep="_")
H3K4me3_gpos100 <- read.table ('H3K4me3_GSM733682cytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
H3K4me3_gpos100$V8 <- paste(H3K4me3_gpos100$V4, H3K4me3_gpos100$V5, H3K4me3_gpos100$V1, sep="_")

CTCF_gpos25 <- read.table ('CTCFcytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
CTCF_gpos25$V8 <- paste(CTCF_gpos25$V4, CTCF_gpos25$V5, CTCF_gpos25$V1, sep="_")
CTCF_gpos50 <- read.table ('CTCFcytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
CTCF_gpos50$V8 <- paste(CTCF_gpos50$V4, CTCF_gpos50$V5, CTCF_gpos50$V1, sep="_")
CTCF_gpos75 <- read.table ('CTCFcytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
CTCF_gpos75$V8 <- paste(CTCF_gpos75$V4, CTCF_gpos75$V5, CTCF_gpos75$V1, sep="_")
CTCF_gpos100 <- read.table ('CTCFcytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
CTCF_gpos100$V8 <- paste(CTCF_gpos100$V4, CTCF_gpos100$V5, CTCF_gpos100$V1, sep="_")

RNAPOLII_gpos25 <- read.table ('RNAPOLIIcytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
RNAPOLII_gpos25$V8 <- paste(RNAPOLII_gpos25$V4, RNAPOLII_gpos25$V5, RNAPOLII_gpos25$V1, sep="_")
RNAPOLII_gpos50 <- read.table ('RNAPOLIIcytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
RNAPOLII_gpos50$V8 <- paste(RNAPOLII_gpos50$V4, RNAPOLII_gpos50$V5, RNAPOLII_gpos50$V1, sep="_")
RNAPOLII_gpos75 <- read.table ('RNAPOLIIcytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
RNAPOLII_gpos75$V8 <- paste(RNAPOLII_gpos75$V4, RNAPOLII_gpos75$V5, RNAPOLII_gpos75$V1, sep="_")
RNAPOLII_gpos100 <- read.table ('RNAPOLIIcytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
RNAPOLII_gpos100$V8 <- paste(RNAPOLII_gpos100$V4, RNAPOLII_gpos100$V5, RNAPOLII_gpos100$V1, sep="_")


H3K9me3_gpos25 <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
H3K9me3_gpos25$V8 <- paste(H3K9me3_gpos25$V4, H3K9me3_gpos25$V5, H3K9me3_gpos25$V1, sep="_")
H3K9me3_gpos50 <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
H3K9me3_gpos50$V8 <- paste(H3K9me3_gpos50$V4, H3K9me3_gpos50$V5, H3K9me3_gpos50$V1, sep="_")
H3K9me3_gpos75 <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
H3K9me3_gpos75$V8 <- paste(H3K9me3_gpos75$V4, H3K9me3_gpos75$V5, H3K9me3_gpos75$V1, sep="_")
H3K9me3_gpos100 <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
H3K9me3_gpos100$V8 <- paste(H3K9me3_gpos100$V4, H3K9me3_gpos100$V5, H3K9me3_gpos100$V1, sep="_")

keep <- rownames(CHD1_gpos25_FDR)
H3K9me3_gpos25 <- H3K9me3_gpos25[keep,]
keep_50 <- rownames(CHD1_gpos50_FDR)
H3K9me3_gpos50 <- H3K9me3_gpos50[keep_50,]
keep_75 <- rownames(CHD1_gpos75_FDR)
H3K9me3_gpos75 <- H3K9me3_gpos75[keep_75,]
keep_100 <- rownames(CHD1_gpos100_FDR)
H3K9me3_gpos100 <- H3K9me3_gpos100[keep_100,]

H3K36me3_gpos25 <- read.table ('H3K36me3_GSM733711cytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
H3K36me3_gpos25$V8 <- paste(H3K36me3_gpos25$V4, H3K36me3_gpos25$V5, H3K36me3_gpos25$V1, sep="_")
H3K36me3_gpos50 <- read.table ('H3K36me3_GSM733711cytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
H3K36me3_gpos50$V8 <- paste(H3K36me3_gpos50$V4, H3K36me3_gpos50$V5, H3K36me3_gpos50$V1, sep="_")
H3K36me3_gpos75 <- read.table ('H3K36me3_GSM733711cytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
H3K36me3_gpos75$V8 <- paste(H3K36me3_gpos75$V4, H3K36me3_gpos75$V5, H3K36me3_gpos75$V1, sep="_")
H3K36me3_gpos100 <- read.table ('H3K36me3_GSM733711cytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
H3K36me3_gpos100$V8 <- paste(H3K36me3_gpos100$V4, H3K36me3_gpos100$V5, H3K36me3_gpos100$V1, sep="_")


H3K4me1_gpos25 <- read.table ('H3K4me1_GSM798322cytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
H3K4me1_gpos25$V8 <- paste(H3K4me1_gpos25$V4, H3K4me1_gpos25$V5, H3K4me1_gpos25$V1, sep="_")
H3K4me1_gpos50 <- read.table ('H3K4me1_GSM798322cytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
H3K4me1_gpos50$V8 <- paste(H3K4me1_gpos50$V4, H3K4me1_gpos50$V5, H3K4me1_gpos50$V1, sep="_")
H3K4me1_gpos75 <- read.table ('H3K4me1_GSM798322cytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
H3K4me1_gpos75$V8 <- paste(H3K4me1_gpos75$V4, H3K4me1_gpos75$V5, H3K4me1_gpos75$V1, sep="_")
H3K4me1_gpos100 <- read.table ('H3K4me1_GSM798322cytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
H3K4me1_gpos100$V8 <- paste(H3K4me1_gpos100$V4, H3K4me1_gpos100$V5, H3K4me1_gpos100$V1, sep="_")

H3K9ac_gpos25 <- read.table ('H3K9ac_GSM733756cytoBand_gpos25.bed_norm_cytob_length', header=FALSE)
H3K9ac_gpos25$V8 <- paste(H3K9ac_gpos25$V4, H3K9ac_gpos25$V5, H3K9ac_gpos25$V1, sep="_")
H3K9ac_gpos50 <- read.table ('H3K9ac_GSM733756cytoBand_gpos50.bed_norm_cytob_length', header=FALSE)
H3K9ac_gpos50$V8 <- paste(H3K9ac_gpos50$V4, H3K9ac_gpos50$V5, H3K9ac_gpos50$V1, sep="_")
H3K9ac_gpos75 <- read.table ('H3K9ac_GSM733756cytoBand_gpos75.bed_norm_cytob_length', header=FALSE)
H3K9ac_gpos75$V8 <- paste(H3K9ac_gpos75$V4, H3K9ac_gpos75$V5, H3K9ac_gpos75$V1, sep="_")
H3K9ac_gpos100 <- read.table ('H3K9ac_GSM733756cytoBand_gpos100.bed_norm_cytob_length', header=FALSE)
H3K9ac_gpos100$V8 <- paste(H3K9ac_gpos100$V4, H3K9ac_gpos100$V5, H3K9ac_gpos100$V1, sep="_")


########################################################################################################
# Generation of an unique table for all the data I read 

total_gpos25_NS <- data.frame(CTCF_gpos25$V7, H3K4me3_gpos25$V7, H3K27ac_gpos25$V7, H3K9ac_gpos25$V7, RNAPOLII_gpos25$V7,
                              CHD1_gpos25_FDR$V7, H3K4me1_gpos25$V7, P300_gpos25$V7, H3K36me3_gpos25$V7,
                              H3K27me3_gpos25$V7, EZH2_gpos25$V7, H3K9me3_gpos25$V7)
dim(total_gpos25_NS)

total_gpos50_NS <- data.frame(CTCF_gpos50$V7, H3K4me3_gpos50$V7, H3K27ac_gpos50$V7, H3K9ac_gpos50$V7, RNAPOLII_gpos50$V7,
                              CHD1_gpos50_FDR$V7, H3K4me1_gpos50$V7, P300_gpos50$V7, H3K36me3_gpos50$V7,
                              H3K27me3_gpos50$V7, EZH2_gpos50$V7, H3K9me3_gpos50$V7)
dim(total_gpos50_NS)

total_gpos75_NS <- data.frame(CTCF_gpos75$V7, H3K4me3_gpos75$V7, H3K27ac_gpos75$V7, H3K9ac_gpos75$V7, RNAPOLII_gpos75$V7,
                              CHD1_gpos75_FDR$V7, H3K4me1_gpos75$V7, P300_gpos75$V7, H3K36me3_gpos75$V7,
                              H3K27me3_gpos75$V7, EZH2_gpos75$V7, H3K9me3_gpos75$V7)
dim(total_gpos75_NS)

total_gpos100_NS <- data.frame(CTCF_gpos100$V7, H3K4me3_gpos100$V7, H3K27ac_gpos100$V7, H3K9ac_gpos100$V7, RNAPOLII_gpos100$V7,
                               CHD1_gpos100_FDR$V7, H3K4me1_gpos100$V7, P300_gpos100$V7, H3K36me3_gpos100$V7,
                               H3K27me3_gpos100$V7, EZH2_gpos100$V7, H3K9me3_gpos100$V7)

dim(total_gpos100_NS)

rownames(total_gpos25_NS) <- CHD1_gpos25_FDR$V8
rownames(total_gpos50_NS) <- CHD1_gpos50_FDR$V8
rownames(total_gpos75_NS) <- CHD1_gpos75_FDR$V8
rownames(total_gpos100_NS) <- CHD1_gpos100_FDR$V8

colnames(total_gpos25_NS) <- c("CTCF", "H3K4me3", "H3K27ac", "H3K9ac", "RNAPOLII", "CHD1", "H3K4me1", "P300",
                               "H3K36me3", "H3K27me3", "EZH2", "H3K9me3")
colnames(total_gpos50_NS) <- c("CTCF", "H3K4me3", "H3K27ac", "H3K9ac", "RNAPOLII", "CHD1", "H3K4me1", "P300",
                               "H3K36me3", "H3K27me3", "EZH2", "H3K9me3")
colnames(total_gpos75_NS) <- c("CTCF", "H3K4me3", "H3K27ac", "H3K9ac", "RNAPOLII", "CHD1", "H3K4me1", "P300",
                               "H3K36me3", "H3K27me3", "EZH2", "H3K9me3")
colnames(total_gpos100_NS) <- c("CTCF", "H3K4me3", "H3K27ac", "H3K9ac", "RNAPOLII", "CHD1", "H3K4me1", "P300",
                                "H3K36me3", "H3K27me3", "EZH2", "H3K9me3")
head(total_gpos25_NS)
colMeans(total_gpos25_NS)  # faster version of apply(scaled.dat, 2, mean) around 0
apply(total_gpos25_NS, 2, sd) #1

total_25_50_75_100_NS <- rbind(total_gpos25_NS, total_gpos50_NS, total_gpos75_NS, total_gpos100_NS)
head(total_25_50_75_100_NS)

#save(total_25_50_75_100_NS, file = "~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/final_objects/Norm_by_cytoband_Length/total_25_50_75_100_NS.RData")
