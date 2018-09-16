# Boxplot showing the distribution of cytogenetic bands length

setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/cytoband_length/")
#
gpos25 <- read.table ("cytoBand_gpos25.bed.length", header=FALSE)
gpos50 <- read.table ("cytoBand_gpos50.bed.length", header=FALSE)
gpos75 <- read.table ("cytoBand_gpos75.bed.length", header=FALSE)
gpos100 <- read.table ("cytoBand_gpos100.bed.length", header=FALSE)

acen <- read.table ("cytoBand_acen.bed.length", header=FALSE)
gneg <- read.table ("cytoBand_gneg.bed.length", header=FALSE)
gvar <- read.table ("cytoBand_gvar.bed.length", header=FALSE)
stalk <- read.table ("cytoBand_stalk.bed.length", header=FALSE)

#head(gpos25)


boxplot (gpos25$V6, gpos50$V6, gpos75$V6, gpos100$V6, acen$V6, gneg$V6, gvar$V6, stalk$V6,
         outline=F, las=2, main="Cytband length", 
         names=c("gpos25", "gpos50", "gpos75", "gpos100", "acen", "gneg", "gvar", "stalk"),
         col = c("gray87", "gray48", "gray34", "gray20","burlywood3", "burlywood3", "burlywood3", "burlywood3"))
