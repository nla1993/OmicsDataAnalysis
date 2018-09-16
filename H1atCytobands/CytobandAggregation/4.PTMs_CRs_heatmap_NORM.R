# Heatmaps of PTMs and CRs
# scaled by cytoband length

library(pheatmap)
library(RColorBrewer)

setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/normalized_by_cytoband_length/")
load("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/final_objects/total_25_50_75_100_NS.RData")

# SCALED DATA #
head(total_25_50_75_100_NS)
PTMs_CRs_NORM <- total_25_50_75_100_NS

#save(PTMs_CRs_NORM, file="~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/final_objects/Norm_by_cytoband_Length/PTMs_CRs_NORM.RData")

PTMs_CRs_scaled <- as.data.frame(apply(total_25_50_75_100_NS, 2, scale))
head(PTMs_CRs_scaled)
rownames(PTMs_CRs_scaled) <- rownames(total_25_50_75_100_NS)
range(PTMs_CRs_scaled)
apply(PTMs_CRs_scaled, 2, range)
hist(as.matrix(PTMs_CRs_scaled))
breaksList = seq(-2.5, 7, by= 0.1)

annotation_prueba <- data.frame(label = c(rep("gpos25", 86), rep("gpos50", 121), rep("gpos75", 89), rep("gpos100", 81)))
head(annotation_prueba)
rownames(annotation_prueba) <- rownames(PTMs_CRs_scaled) # check out the row names of annotation

# 1. Manhattan clustering distance in rows
# pheatmap(PTMs_CRs_scaled, cluster_rows = T, cluster_cols = T, scale = "none", 
#          breaks = breaksList, fontsize_row = 3, show_rownames = F, clustering_distance_rows = "manhattan",
#          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
#          annotation_row = annotation_prueba, cutree_cols = 6)

# 2. Correlation clustering distance in rows

hm1 <- pheatmap(PTMs_CRs_scaled, cluster_rows = T, cluster_cols = T, scale = "none", 
                breaks = breaksList, fontsize_row = 3, show_rownames = F, clustering_distance_rows = "correlation",
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                annotation_row = annotation_prueba, cutree_cols = 6, cutree_rows = 8)

# 3. Euclidean clustering distance in rows
# pheatmap(PTMs_CRs_scaled, cluster_rows = T, cluster_cols = T, scale = "none", 
#          breaks = breaksList, fontsize_row = 3, show_rownames = F, clustering_distance_rows = "euclidean",
#          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
#          annotation_row = annotation_prueba, cutree_cols = 6)

# 4. Minkowski clustering distance in rows
# pheatmap(PTMs_CRs_scaled, cluster_rows = T, cluster_cols = T, scale = "none", 
#          breaks = breaksList, fontsize_row = 3, show_rownames = F, clustering_distance_rows = "minkowski",
#          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
#          annotation_row = annotation_prueba, cutree_cols = 6)

#dev.off()