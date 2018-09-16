# Boxplot of gene expression groups

setwd("/Users/Natalia/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/HeLa_expression_RNA-Seq/")
Not <- read.table ('HeLaS3_RNAseq_sorted_NotExpressed.strand.bed', header=FALSE)
xa <- read.table ('xaa_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xb <- read.table ('xab_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xc <- read.table ('xac_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xd <- read.table ('xad_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xe <- read.table ('xae_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xf <- read.table ('xaf_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xg <- read.table ('xag_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xh <- read.table ('xah_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xi <- read.table ('xai_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)
xj <- read.table ('xaj_HeLaS3_RNAseq_sorted_Expressed.strand.bed', header=FALSE)

#pdf("Boxplot_Expression_Groups_HeLaRNAseq.pdf")
boxplot (Not$V6, xa$V6, xb$V6, xc$V6, xd$V6, xe$V6, xf$V6, xg$V6, xh$V6, xi$V6, xj$V6,  outline=FALSE, las=2, 
         ylab="Expression (RPKM)", names=c("Not", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
         main= "Gene Expression Groups Distribution")
#dev.off()

