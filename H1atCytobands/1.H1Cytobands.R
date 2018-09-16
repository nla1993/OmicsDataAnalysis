## 1. GC content at cytobands (without Ns)
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/cytoBand_HeLa/andrea_cytoband/GCcontent/")

acen <- read.table ('cytoBand_acen.GCcontent_NotNs', header=FALSE)
gneg <- read.table ('cytoBand_gneg.GCcontent_NotNs', header=FALSE)
gpos100 <- read.table ('cytoBand_gpos100.GCcontent_NotNs', header=FALSE)
gpos75 <- read.table ('cytoBand_gpos75.GCcontent_NotNs', header=FALSE)
gpos50 <- read.table ('cytoBand_gpos50.GCcontent_NotNs', header=FALSE)
gpos25 <- read.table ('cytoBand_gpos25.GCcontent_NotNs', header=FALSE)
gvar <- read.table ('cytoBand_gvar.GCcontent_NotNs', header=FALSE)
stalk <- read.table ('cytoBand_stalk.GCcontent_NotNs', header=FALSE)

boxplot (gneg$V4*100, gpos25$V4*100, gpos50$V4*100, gpos75$V4*100, gpos100$V4*100, acen$V4*100, gvar$V4*100, 
         outline=FALSE, las=2, ylab="GC content (%)", main="GC content \nat cytoband groups", ylim = c(30,60),
         names=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"), ylim=c(0,100),
         col = c("burlywood3", "gray87", "gray48", "gray34", "gray20", "burlywood3", "burlywood3"))

# 2. PTMs abundance at cytobands
## H3K9me3 
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/cytoBand_HeLa/andrea_cytoband/H3K9me3/")

## broadPeaks
acen <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_acen.bed', header=FALSE)
gneg <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gneg.bed', header=FALSE)
gpos100 <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gpos100.bed', header=FALSE)
gpos75 <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gpos75.bed', header=FALSE)
gpos50 <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gpos50.bed', header=FALSE)
gpos25 <- read.table ('H3K9me3_HeLa_broadPeaksVscytoBand_gpos25.bed', header=FALSE)
gvar <- read.table('H3K9me3_HeLa_broadPeaksVscytoBand_gvar.bed', header=FALSE)

boxplot (gneg$V6, gpos25$V6, gpos50$V6, gpos75$V6, gpos100$V6, acen$V6, gvar$V6,
         outline=FALSE, las=2, main="H3K9me3 abundance \nat cytoband groups", 
         names=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"),
         col = c("burlywood3", "gray87", "gray48", "gray34", "gray20", "burlywood3", "burlywood3"))

## H3K27me3
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/cytoBand_HeLa/H3K27me3/H3K27me3_GSM733696/")
gneg_H3K27me3 <- read.table ('H3K27me3_GSM733696cytoBand_gneg.bed', header=FALSE)
gpos25_H3K27me3 <- read.table ('H3K27me3_GSM733696cytoBand_gpos25.bed', header=FALSE)
gpos50_H3K27me3 <- read.table ('H3K27me3_GSM733696cytoBand_gpos50.bed', header=FALSE)
gpos75_H3K27me3 <- read.table ('H3K27me3_GSM733696cytoBand_gpos75.bed', header=FALSE)
gpos100_H3K27me3 <- read.table ('H3K27me3_GSM733696cytoBand_gpos100.bed', header=FALSE)
gvar_H3K27me3 <- read.table ('H3K27me3_GSM733696cytoBand_gvar.bed', header=FALSE)
acen_H3K27me3 <- read.table ('H3K27me3_GSM733696cytoBand_acen.bed', header=FALSE)
stalk_H3K27me3 <- read.table ('H3K27me3_GSM733696cytoBand_stalk.bed', header=FALSE)

sapply(gneg_H3K27me3, class)

head(gneg_H3K27me3)
gneg_H3K27me3 <- gneg_H3K27me3[, "V6"]
gpos25_H3K27me3 <- gpos25_H3K27me3[, "V6"]
gpos50_H3K27me3 <- gpos50_H3K27me3[, "V6"]
gpos75_H3K27me3 <- gpos75_H3K27me3[, "V6"]
gpos100_H3K27me3 <- gpos100_H3K27me3[, "V6"]
gvar_H3K27me3 <- gvar_H3K27me3[, "V6"]
acen_H3K27me3 <- acen_H3K27me3[, "V6"]
stalk_H3K27me3 <- stalk_H3K27me3[, "V6"]

boxplot(gneg_H3K27me3, gpos25_H3K27me3, gpos50_H3K27me3, gpos75_H3K27me3, gpos100_H3K27me3,
        acen_H3K27me3, gvar_H3K27me3, outline = FALSE, las=2, main = "H3K27me3 abundance \nat cytoband groups",
        names=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"), 
        col = c("burlywood3", "gray87", "gray48", "gray34", "gray20", "burlywood3", "burlywood3")) 

## H3K4me3
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/cytoBand_HeLa/H3K4me3/H3K4me3_GSM733682/")
gneg_H3K4me3 <- read.table ('H3K4me3_GSM733682cytoBand_gneg.bed', header=FALSE)
gpos25_H3K4me3 <- read.table ('H3K4me3_GSM733682cytoBand_gpos25.bed', header=FALSE)
gpos50_H3K4me3 <- read.table ('H3K4me3_GSM733682cytoBand_gpos50.bed', header=FALSE)
gpos75_H3K4me3 <- read.table ('H3K4me3_GSM733682cytoBand_gpos75.bed', header=FALSE)
gpos100_H3K4me3 <- read.table ('H3K4me3_GSM733682cytoBand_gpos100.bed', header=FALSE)
gvar_H3K4me3 <- read.table ('H3K4me3_GSM733682cytoBand_gvar.bed', header=FALSE)
acen_H3K4me3 <- read.table ('H3K4me3_GSM733682cytoBand_acen.bed', header=FALSE)
stalk_H3K4me3 <- read.table ('H3K4me3_GSM733682cytoBand_stalk.bed', header=FALSE)

sapply(gneg_H3K4me3, class)
sapply(gpos25_H3K4me3, class)
sapply(gpos50_H3K4me3, class)
sapply(gpos75_H3K4me3, class)
sapply(gpos100_H3K4me3, class)
sapply(gvar_H3K4me3, class)
sapply(acen_H3K4me3, class)
sapply(stalk_H3K4me3, class)


head(gneg_H3K4me3)
gneg_H3K4me3 <- gneg_H3K4me3[, "V6"]
gpos25_H3K4me3 <- gpos25_H3K4me3[, "V6"]
gpos50_H3K4me3 <- gpos50_H3K4me3[, "V6"]
gpos75_H3K4me3 <- gpos75_H3K4me3[, "V6"]
gpos100_H3K4me3 <- gpos100_H3K4me3[, "V6"]
gvar_H3K4me3 <- gvar_H3K4me3[, "V6"]
acen_H3K4me3 <- acen_H3K4me3[, "V6"]
stalk_H3K4me3 <- stalk_H3K4me3[, "V6"]

boxplot(gneg_H3K4me3, gpos25_H3K4me3, gpos50_H3K4me3, gpos75_H3K4me3, gpos100_H3K4me3,
        acen_H3K4me3, gvar_H3K4me3, outline = FALSE, las=2, main = "H3K4me3 abundance \nat cytoband groups",
        names=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"),
        col = c("burlywood3", "gray87", "gray48", "gray34", "gray20", "burlywood3", "burlywood3")) 


# 3. H1 variants abundance at cytobands
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/cytoBand_HeLa/andrea_cytoband/vsH1/")

H1.2_acen <- read.table ('H1_2_InputMergedSubtracted_fromBam.bdgVSacen', header=FALSE)
H1X_acen <- read.table ('H1X_InputMergedSubtracted_fromBam.bdgVSacen', header=FALSE)
HA_0_acen <- read.table ('0_HA_InputMergedSubtracted_fromBam.bdgVSacen', header=FALSE)
HA_2_acen <- read.table ('2_HA_InputMergedSubtracted_fromBam.bdgVSacen', header=FALSE)

H1.2_stalk <- read.table ('H1_2_InputMergedSubtracted_fromBam.bdgVSstalk', header=FALSE)
H1X_stalk <- read.table ('H1X_InputMergedSubtracted_fromBam.bdgVSstalk', header=FALSE)
HA_0_stalk <- read.table ('0_HA_InputMergedSubtracted_fromBam.bdgVSstalk', header=FALSE)
HA_2_stalk <- read.table ('2_HA_InputMergedSubtracted_fromBam.bdgVSstalk', header=FALSE)

H1.2_gvar <- read.table ('H1_2_InputMergedSubtracted_fromBam.bdgVSgvar', header=FALSE)
H1X_gvar <- read.table ('H1X_InputMergedSubtracted_fromBam.bdgVSgvar', header=FALSE)
HA_0_gvar <- read.table ('0_HA_InputMergedSubtracted_fromBam.bdgVSgvar', header=FALSE)
HA_2_gvar <- read.table ('2_HA_InputMergedSubtracted_fromBam.bdgVSgvar', header=FALSE)

H1.2_gneg <- read.table ('H1_2_InputMergedSubtracted_fromBam.bdgVSgneg', header=FALSE)
H1X_gneg <- read.table ('H1X_InputMergedSubtracted_fromBam.bdgVSgneg', header=FALSE)
HA_0_gneg <- read.table ('0_HA_InputMergedSubtracted_fromBam.bdgVSgneg', header=FALSE)
HA_2_gneg <- read.table ('2_HA_InputMergedSubtracted_fromBam.bdgVSgneg', header=FALSE)

H1.2_gpos25 <- read.table ('H1_2_InputMergedSubtracted_fromBam.bdgVSgpos25', header=FALSE)
H1X_gpos25 <- read.table ('H1X_InputMergedSubtracted_fromBam.bdgVSgpos25', header=FALSE)
HA_0_gpos25 <- read.table ('0_HA_InputMergedSubtracted_fromBam.bdgVSgpos25', header=FALSE)
HA_2_gpos25 <- read.table ('2_HA_InputMergedSubtracted_fromBam.bdgVSgpos25', header=FALSE)

H1.2_gpos25$V6 <- as.numeric(as.character(H1.2_gpos25$V6))
H1X_gpos25$V6 <- as.numeric(as.character(H1X_gpos25$V6))
HA_0_gpos25$V6 <- as.numeric(as.character(HA_0_gpos25$V6))
HA_2_gpos25$V6 <- as.numeric(as.character(HA_2_gpos25$V6))

H1.2_gpos50 <- read.table ('H1_2_InputMergedSubtracted_fromBam.bdgVSgpos50', header=FALSE)
H1X_gpos50 <- read.table ('H1X_InputMergedSubtracted_fromBam.bdgVSgpos50', header=FALSE)
HA_0_gpos50 <- read.table ('0_HA_InputMergedSubtracted_fromBam.bdgVSgpos50', header=FALSE)
HA_2_gpos50 <- read.table ('2_HA_InputMergedSubtracted_fromBam.bdgVSgpos50', header=FALSE)

H1.2_gpos75 <- read.table ('H1_2_InputMergedSubtracted_fromBam.bdgVSgpos75', header=FALSE)
H1X_gpos75 <- read.table ('H1X_InputMergedSubtracted_fromBam.bdgVSgpos75', header=FALSE)
HA_0_gpos75 <- read.table ('0_HA_InputMergedSubtracted_fromBam.bdgVSgpos75', header=FALSE)
HA_2_gpos75 <- read.table ('2_HA_InputMergedSubtracted_fromBam.bdgVSgpos75', header=FALSE)

H1.2_gpos100 <- read.table ('H1_2_InputMergedSubtracted_fromBam.bdgVSgpos100', header=FALSE)
H1X_gpos100 <- read.table ('H1X_InputMergedSubtracted_fromBam.bdgVSgpos100', header=FALSE)
HA_0_gpos100 <- read.table ('0_HA_InputMergedSubtracted_fromBam.bdgVSgpos100', header=FALSE)
HA_2_gpos100 <- read.table ('2_HA_InputMergedSubtracted_fromBam.bdgVSgpos100', header=FALSE)

H1X_gpos25$V6 <- as.numeric(as.character(H1X_gpos25$V6))

boxplot (H1.2_gneg$V6, H1.2_gpos25$V6, H1.2_gpos50$V6, H1.2_gpos75$V6, H1.2_gpos100$V6, H1.2_acen$V6, H1.2_gvar$V6,
         outline=FALSE, las=2, main= "H1.2-endo abundance \nat cytoband groups",
         names=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"),
         col = c("burlywood3", "gray87", "gray48", "gray34", "gray20", "burlywood3", "burlywood3"))
boxplot (H1X_gneg$V6, H1X_gpos25$V6, H1X_gpos50$V6, H1X_gpos75$V6, H1X_gpos100$V6, H1X_acen$V6, H1X_gvar$V6,
         outline=FALSE, las=2, main= "H1.X-endo abundance \nat cytoband groups",
         names=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"),
         col = c("burlywood3", "gray87", "gray48", "gray34", "gray20", "burlywood3", "burlywood3"))
boxplot (HA_2_gneg$V6, HA_2_gpos25$V6, HA_2_gpos50$V6, HA_2_gpos75$V6, HA_2_gpos100$V6, HA_2_acen$V6, HA_2_gvar$V6,
         outline=FALSE, las=2 , main= "H1.2-HA abundance \nat cytoband groups",
         names=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"),
         col = c("burlywood3", "gray87", "gray48", "gray34", "gray20", "burlywood3", "burlywood3"))
boxplot (HA_0_gneg$V6, HA_0_gpos25$V6, HA_0_gpos50$V6, HA_0_gpos75$V6, HA_0_gpos100$V6, HA_0_acen$V6, HA_0_gvar$V6,
         outline=FALSE, las=2, main= "H1.0-HA abundance \nat cytoband groups",
         names=c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"),
         col = c("burlywood3", "gray87", "gray48", "gray34", "gray20", "burlywood3", "burlywood3"))
