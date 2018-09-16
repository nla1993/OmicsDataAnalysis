# LADs
# Boxplot showing the occupancy of H1 variants at LADS, pericentromeric regions, accesible DNA
# regions (FAIRE-seq), and enhancers.
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/180420_LADS/")

H1_2_LADS_mean <- read.table ('H1_2_LADS_mean', header=FALSE)
H1X_LADS_mean <- read.table ('H1X_LADS_mean', header=FALSE)
H1.0_HA_LADS_mean <- read.table ('0_HA_LADS_mean', header=FALSE)
H1.2_HA_LADS_mean <- read.table ('2_HA_LADS_mean', header=FALSE) 


final_table <- cbind(H1_2_LADS_mean, H1X_LADS_mean, H1.0_HA_LADS_mean, H1.2_HA_LADS_mean)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "LADs", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                    "sandybrown"))
abline(h = 0, col = "grey", lwd = 1) 

# Centromeres
# Boxplot showing the occupancy of H1 variants at centromeric region in p&q arms
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/centromeres/H1_centromere_abundance_heatmap/")

# p arm data
HA_0_p <- read.table("0_HA_InputMergedSubtracted_fromBam_abundance_p_centrom", header = FALSE, sep = "\t")
HA_2_p <- read.table("2_HA_InputMergedSubtracted_fromBam_abundance_p_centrom", header = FALSE, sep = "\t")
H1_X_p <- read.table("H1X_InputMergedSubtracted_fromBam_abundance_p_centrom", header = FALSE, sep = "\t")
H1_2_p <- read.table("H1_2_InputMergedSubtracted_fromBam_abundance_p_centrom", header = FALSE, sep = "\t")

colnames(HA_0_p) <- c("V1", "HA_0")
colnames(HA_2_p) <- c("V1", "HA_2")
colnames(H1_X_p) <- c("V1", "H1_X")
colnames(H1_2_p) <- c("V1", "H1_2")

data <- merge(HA_0_p, HA_2_p, by="V1")
data <- merge(data, H1_X_p, by="V1")
data <- merge(data, H1_2_p, by="V1")

head(data)

rownames(data) <- data$V1
data_p <- data[, -1]
head(data_p)
dim(data_p)

# q arm data
HA_0_q <- read.table("0_HA_InputMergedSubtracted_fromBam_abundance_q_centrom", header = FALSE, sep = "\t")
HA_2_q <- read.table("2_HA_InputMergedSubtracted_fromBam_abundance_q_centrom", header = FALSE, sep = "\t")
H1_X_q <- read.table("H1X_InputMergedSubtracted_fromBam_abundance_q_centrom", header = FALSE, sep = "\t")
H1_2_q <- read.table("H1_2_InputMergedSubtracted_fromBam_abundance_q_centrom", header = FALSE, sep = "\t")

colnames(HA_0_q) <- c("V1", "HA_0")
colnames(HA_2_q) <- c("V1", "HA_2")
colnames(H1_X_q) <- c("V1", "H1_X")
colnames(H1_2_q) <- c("V1", "H1_2")

data <- merge(HA_0_q, HA_2_q, by="V1")
data <- merge(data, H1_X_q, by="V1")
data <- merge(data, H1_2_q, by="V1")

head(data)

rownames(data) <- data$V1
data_q <- data[, -1]
head(data_q)
dim(data_q)

data_p_ordered <- data_p[, c("H1_2", "H1_X", "HA_0", "HA_2")]
colnames(data_p_ordered) <- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")
boxplot(data_p_ordered, outline = FALSE, las=2, main = "centromeric position \nat p arms", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                                                 "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)

data_q_ordered <- data_q[, c("H1_2", "H1_X", "HA_0", "HA_2")]
colnames(data_q_ordered) <- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(data_q_ordered, outline = FALSE, las=2, main = "centromeric position \nat q arms", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                                                 "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)

# Abundance of H1 variants at accesible DNA regions (FAIRE seq data)
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/FAIRE-Seq_HeLa_S3/")

H1_2_mean <- read.table ('H1_2_FAIRE_peaks_overlap_mean', header=FALSE)
H1X_mean <- read.table ('H1X_FAIRE_peaks_overlap_mean', header=FALSE)

H1_0_HA_mean <- read.table ('0_HA_FAIRE_peaks_overlap_mean', header=FALSE)
H1_2_HA_mean <- read.table ('2_HA_FAIRE_peaks_overlap_mean', header=FALSE)


final_table <- cbind(H1_2_mean, H1X_mean, H1_0_HA_mean, H1_2_HA_mean)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")
boxplot(final_table, outline = FALSE, las=2, main = "Accessible DNA regions \n (FAIRE peaks)", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                                                      "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)

# Abundance of H1 variants at enhancers
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/180319_chIP_seq/enhancers_HeLa-S3")

H1_2_enhancers_overlap<- read.table ('H1_2_enhancers_peaks_overlap', header=FALSE)
H1X_enhancers_overlap<- read.table ('H1X_enhancers_peaks_overlap', header=FALSE)

H1.0_HA_enhancers_overlap<- read.table ('0_HA_enhancers_peaks_overlap', header=FALSE)
H1.2_HA_enhancers_overlap<- read.table ('2_HA_enhancers_peaks_overlap', header=FALSE)


final_table <- cbind (H1_2_enhancers_overlap$V4, H1X_enhancers_overlap$V4, H1.0_HA_enhancers_overlap$V4,
                      H1.2_HA_enhancers_overlap$V4)
colnames(final_table)<- c("H1.2-endo", "H1.X-endo", "H1.0-HA", "H1.2-HA")

boxplot(final_table, outline = FALSE, las=2, main = "Enhancers", col = c("sandybrown", "lightblue1", "lightpink", 
                                                                         "sandybrown"))
abline(h = 0, col = "grey", lwd = 1)

