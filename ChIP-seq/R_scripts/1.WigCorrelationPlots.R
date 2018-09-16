# Wig Correlation Plots
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/Wig_correlation")

# 1. Heatmap of Wig files previous to Input subtraction
all_data <- read.table("180406_wigCorrelation_table.csv", header = TRUE, sep = ";")
colname_all_data <- c("wt-input", "HA-input", "fHA-input", "H1.0-HA", "H1.0-flagHA", "H1.2-endo", "H1.2-HA", "H1.2-flagHA", "H1.4-flagHA",
                      "H1.X-endo", "H1.X-flagHA")
rownames(all_data) <- colname_all_data
all_mydata <- all_data[,-1]
colnames(all_mydata) <- colname_all_data
head(all_mydata)

heatmap(as.matrix(all_mydata),cexRow=0.85, cexCol=0.85, scale="none",
        col=paste("grey",99:0,sep=""))

# 2. Heatmap of Wig correlation after Input merged subtraction
mydata <- read.table("180411_WigCorrelation_InputSubstracted_fromBam.csv", header = TRUE, sep = ";")
colname <- c("H1.0-HA", "H1.0-flagHA", "H1.2-endo", "H1.2-HA", "H1.2-flagHA", "H1.4-flagHA",
             "H1.X-endo", "H1.X-flagHA")

rownames(mydata) <- colname
mydata <- mydata[,-1]
colnames(mydata) <- colname
head(mydata)

heatmap(as.matrix(mydata),cexRow=0.85, cexCol=0.85, scale="none",
        col=paste("grey",99:0,sep=""))

# 3. Representation of Correlation values (Input merged subtraction)
cormat <- round(mydata,2)
head(cormat)
class(cormat)
cormat <- as.matrix(cormat)

library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
sapply(melted_cormat, class)

# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.6, limit = c(0.4,1), space = "Lab", 
                       name="Wig\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4)

# 4. Samples we will explore in further analyses
mydata_2 <- mydata[c(1,3,4,7), c(1,3,4,7)]

heatmap(as.matrix(mydata_2),cexRow=0.85, cexCol=0.85, scale="none",
        col=paste("grey",99:0,sep=""))

dev.off()
