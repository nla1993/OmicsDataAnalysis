# Cytoband aggregation
## 1. Table of H1 abundance at cytobands

### Reading the files
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_hela_aggregation/H1variantsOverlap/")

H1.2_gpos25 <- read.table ('H1_2_cytoBand_gpos25_overlap.bed', header=FALSE)
H1.2_gpos25$V6 <- as.numeric(as.character(H1.2_gpos25$V6))
H1.2_gpos25$V7 <- paste(H1.2_gpos25$V4, H1.2_gpos25$V5, H1.2_gpos25$V1, sep="_")
H1.2_gpos25 <- na.omit(H1.2_gpos25)
H1.2_gpos50 <- read.table ('H1_2_cytoBand_gpos50_overlap.bed', header=FALSE)
H1.2_gpos50$V6 <- as.numeric(as.character(H1.2_gpos50$V6))
H1.2_gpos50$V7 <- paste(H1.2_gpos50$V4, H1.2_gpos50$V5, H1.2_gpos50$V1, sep="_")
H1.2_gpos75 <- read.table ('H1_2_cytoBand_gpos75_overlap.bed', header=FALSE)
H1.2_gpos75$V6 <- as.numeric(as.character(H1.2_gpos75$V6))
H1.2_gpos75$V7 <- paste(H1.2_gpos75$V4, H1.2_gpos75$V5, H1.2_gpos75$V1, sep="_")
H1.2_gpos100 <- read.table ('H1_2_cytoBand_gpos100_overlap.bed', header=FALSE)
H1.2_gpos100$V6 <- as.numeric(as.character(H1.2_gpos100$V6))
H1.2_gpos100$V7 <- paste(H1.2_gpos100$V4, H1.2_gpos100$V5, H1.2_gpos100$V1, sep="_")


H1X_gpos25 <- read.table ('H1X_cytoBand_gpos25_overlap.bed', header=FALSE)
H1X_gpos25 <- na.omit(H1X_gpos25)
H1X_gpos25$V6 <- as.numeric(as.character(H1X_gpos25$V6))
H1X_gpos25$V7 <- paste(H1X_gpos25$V4, H1X_gpos25$V5, H1X_gpos25$V1, sep="_")


keep <- rownames(H1.2_gpos25)
H1X_gpos25 <- H1X_gpos25[keep,]

H1X_gpos50 <- read.table ('H1X_cytoBand_gpos50_overlap.bed', header=FALSE)
H1X_gpos50$V6 <- as.numeric(as.character(H1X_gpos50$V6))
H1X_gpos50$V7 <- paste(H1X_gpos50$V4, H1X_gpos50$V5, H1X_gpos50$V1, sep="_")
H1X_gpos75 <- read.table ('H1X_cytoBand_gpos75_overlap.bed', header=FALSE)
H1X_gpos75$V6 <- as.numeric(as.character(H1X_gpos75$V6))
H1X_gpos75$V7 <- paste(H1X_gpos75$V4, H1X_gpos75$V5, H1X_gpos75$V1, sep="_")
H1X_gpos100 <- read.table ('H1X_cytoBand_gpos100_overlap.bed', header=FALSE)
H1X_gpos100$V6 <- as.numeric(as.character(H1X_gpos100$V6))
H1X_gpos100$V7 <- paste(H1X_gpos100$V4, H1X_gpos100$V5, H1X_gpos100$V1, sep="_")

dim(H1.2_gpos25)


total_gpos25 <- data.frame(H1.2_gpos25$V6, H1X_gpos25$V6)
dim(total_gpos25)

total_gpos50 <- data.frame(H1.2_gpos50$V6, H1X_gpos50$V6)
dim(total_gpos50)

total_gpos75 <- data.frame(H1.2_gpos75$V6, H1X_gpos75$V6)
dim(total_gpos75)

total_gpos100 <- data.frame(H1.2_gpos100$V6, H1X_gpos100$V6)
dim(total_gpos100)

rownames(total_gpos25) <- H1.2_gpos25$V7
rownames(total_gpos50) <- H1.2_gpos50$V7
rownames(total_gpos75) <- H1.2_gpos75$V7
rownames(total_gpos100) <- H1.2_gpos100$V7

colnames(total_gpos25) <- c("H1.2_HeLa", "H1X_HeLa")
colnames(total_gpos50) <- c("H1.2_HeLa", "H1X_HeLa")
colnames(total_gpos75) <- c("H1.2_HeLa", "H1X_HeLa")
colnames(total_gpos100) <- c("H1.2_HeLa", "H1X_HeLa")

head(total_gpos25)
sapply(total_gpos25, class)
# colMeans(total_gpos25)  # faster version of apply(scaled.dat, 2, mean) around 0
# apply(total_gpos25_NS, 2, sd) #1

total_25_50_75_100 <- rbind(total_gpos25, total_gpos50, total_gpos75, total_gpos100)
total_25_50_75_100 <- total_25_50_75_100*1000
#total_25_50_75_100 <- na.omit(total_25_50_75_100)

H1_total_25_50_75_100 <- total_25_50_75_100
head(H1_total_25_50_75_100)

H1variants <- H1_total_25_50_75_100

# scale the data
H1variants_hela_scaled <-  as.data.frame(apply(H1variants, 2, scale))
head(H1variants_hela_scaled)
rownames(H1variants_hela_scaled) <- rownames(H1variants)
head(H1variants_hela_scaled)

apply(H1variants_hela_scaled, 2, range)
table(is.na(H1variants_hela_scaled))

