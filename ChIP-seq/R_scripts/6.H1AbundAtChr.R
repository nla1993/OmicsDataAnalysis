# H1 abundance at chromosome level
# Just for samples for the article, and some chromosomes (not all of them)

setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/180417_h1_abundance_chr_heatmap")
fHA_0 <- read.table("0_fHA_InputMergedSubtracted_fromBam.bdg_norm_val_chr", header = FALSE, sep = "\t")
HA_0 <- read.table("0_HA_InputMergedSubtracted_fromBam.bdg_norm_val_chr", header = FALSE, sep = "\t")
fHA_2 <- read.table("2_fHA_InputMergedSubtracted_fromBam.bdg_norm_val_chr", header = FALSE, sep = "\t")
HA_2 <- read.table("2_HA_InputMergedSubtracted_fromBam.bdg_norm_val_chr", header = FALSE, sep = "\t")
fHA_4 <- read.table("4_fHA_InputMergedSubtracted_fromBam.bdg_norm_val_chr", header = FALSE, sep = "\t")
fHA_X <- read.table("X_fHA_InputMergedSubtracted_fromBam.bdg_norm_val_chr", header = FALSE, sep = "\t")
H1_X <- read.table("H1X_InputMergedSubtracted_fromBam.bdg_norm_val_chr", header = FALSE, sep = "\t")
H1_2 <- read.table("H1_2_InputMergedSubtracted_fromBam.bdg_norm_val_chr", header = FALSE, sep = "\t")

head(fHA_0)
head(HA_0)
head(fHA_2)
head(HA_2)
head(H1_X)

colnames(fHA_0) <- c("V1", "fHA_0")
colnames(HA_0) <- c("V1", "HA_0")
colnames(fHA_2) <- c("V1", "fHA_2")
colnames(HA_2) <- c("V1", "HA_2")
colnames(fHA_4) <- c("V1", "fHA_4")
colnames(fHA_X) <- c("V1", "fHA_X")
colnames(H1_X) <- c("V1", "H1_X")
colnames(H1_2) <- c("V1", "H1_2")

data <- merge(fHA_0, HA_0, by="V1")
data <- merge(data, fHA_2, by="V1")
data <- merge(data, HA_2, by="V1")
data <- merge(data, fHA_4, by="V1")
data <- merge(data, fHA_X, by="V1")
data <- merge(data, H1_X, by="V1")
data <- merge(data, H1_2, by="V1")
head(data)

rownames(data) <- data$V1
data <- data[, -1]
head(data)
dim(data)

data_1 <- data[-23, ]
head(data_1)
dim(data_1)
data_2 <- data_1[-24, ]
head(data_2)
dim(data_2)

head(data_2)
mydata <- data_2
mydata <- mydata * 1000
head(mydata)

colname <- c("H1.0-flagHA", "H1.0-HA", "H1.2-flagHA", "H1.2-HA", "H1.4-flagHA", 
             "H1.X-flagHA", "H1.X-endo", "H1.2-endo")
# 
# rownames(mydata) <- colname
# mydata <- mydata[,-1]

colnames(mydata) <- colname
head(mydata)

#ordering
mydata <- mydata[,c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA", "H1.0-flagHA", "H1.2-flagHA",
                    "H1.4-flagHA", "H1.X-flagHA")]

mydata <- data_2[c(9,15,2,18,1,5,11,19,6,21,23),c("H1_X","H1_2","HA_2", "HA_0")]

head(mydata)

colname <- c("H1.X-endo", "H1.2-endo", "H1.2-HA", "H1.0-HA")

colnames(mydata) <- colname
head(mydata)
mydata <- mydata * 1000

cormat <- round(mydata,2)
head(cormat)
class(cormat)
cormat <- as.matrix(cormat)

library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
sapply(melted_cormat, class)

# library(ggplot2)
# ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile()

# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-2,4.5), space = "Lab", 
                       name="H1 variants\n Abundance (x1000)") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4)
