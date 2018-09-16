###############################################
setwd("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/normalized_by_cytoband_length/")

source("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/final scripts/Norm_by_cytoband_Length/2.PTMs_CRs_heatmap_NORM.R")
#load("~/Documents/OneDrive - Universitat de Vic/TFM/Transciptomics/CytoBand_HeLa_aggregation/final_objects/H1variants_hela_scaled.RData")


# H1.2 and H1X alone but scaled before together with the PTMs & CRs
order <- hm1$tree_row$order
head(order)

# Scaled by their own
# Data from script H1variants_abundance_cytoband_aggregation
head(H1variants_hela_scaled)
H1variants_hela_scaled_order <- H1variants_hela_scaled
head(H1variants_hela_scaled_order)
H1variants_hela_scaled_order$V3 <- rownames(H1variants_hela_scaled_order)
rownames(H1variants_hela_scaled_order) <- NULL
head(H1variants_hela_scaled_order)
H1variants_hela_scaled_order$V4 <- rownames(H1variants_hela_scaled_order)

H1variants_hela_ordered <- H1variants_hela_scaled_order[match(order, H1variants_hela_scaled_order$V4),]
head(H1variants_hela_ordered)

# defining the groups
prueba <- H1variants_hela_ordered
head(prueba)
prueba$V5 <- c(1:length(prueba$V4))

prueba[prueba$V3=="q31.33_gpos75_chr7", ] #group 1
prueba[prueba$V3=="q23.1_gpos75_chr16", ] #group 2
prueba[prueba$V3=="q26.33_gpos75_chr3", ]
prueba[prueba$V3=="q13.2_gpos50_chrX", ] #group 4
prueba[prueba$V3=="q28.3_gpos100_chr4", ]
prueba[prueba$V3=="q23.31_gpos75_chr10", ]

rownames(prueba) <- prueba$V5
head(prueba)
group_1 <- prueba[c(1:69),]
group_2 <- prueba[c(70:108),]
group_3 <- prueba[c(109:131),]
group_4 <- prueba[c(132:175),]
group_5 <- prueba[c(176:226),]
group_6 <- prueba[c(226:243),]
group_7 <- prueba[c(244:312),]
group_8 <- prueba[c(313:377),]


head(group_1)
library(stringr)
dim(group_1)
table(str_count(group_1$V3, "gpos25"))
table(str_count(group_1$V3, "gpos50"))
table(str_count(group_1$V3, "gpos75"))
table(str_count(group_1$V3, "gpos100"))

group_1_table <- data.frame(c(5,13,18,33))

table(str_count(group_2$V3, "gpos25"))
table(str_count(group_2$V3, "gpos50"))
table(str_count(group_2$V3, "gpos75"))
table(str_count(group_2$V3, "gpos100"))

group_2_table <- data.frame(c(13,16,7,3))

table(str_count(group_3$V3, "gpos25"))
table(str_count(group_3$V3, "gpos50"))
table(str_count(group_3$V3, "gpos75"))
table(str_count(group_3$V3, "gpos100"))

group_3_table <- data.frame(c(8,13,2,0))

table(str_count(group_4$V3, "gpos25"))
table(str_count(group_4$V3, "gpos50"))
table(str_count(group_4$V3, "gpos75"))
table(str_count(group_4$V3, "gpos100"))

group_4_table <- data.frame(c(8,16,15,5))

table(str_count(group_5$V3, "gpos25"))
table(str_count(group_5$V3, "gpos50"))
table(str_count(group_5$V3, "gpos75"))
table(str_count(group_5$V3, "gpos100"))

group_5_table <- data.frame(c(32,17,1,1))

table(str_count(group_6$V3, "gpos25"))
table(str_count(group_6$V3, "gpos50"))
table(str_count(group_6$V3, "gpos75"))
table(str_count(group_6$V3, "gpos100"))

group_6_table <- data.frame(c(3,7,1,7))

table(str_count(group_7$V3, "gpos25"))
table(str_count(group_7$V3, "gpos50"))
table(str_count(group_7$V3, "gpos75"))
table(str_count(group_7$V3, "gpos100"))

group_7_table <- data.frame(c(12,17,21,19))

table(str_count(group_8$V3, "gpos25"))
table(str_count(group_8$V3, "gpos50"))
table(str_count(group_8$V3, "gpos75"))
table(str_count(group_8$V3, "gpos100"))

group_8_table <- data.frame(c(6,22,24,13))

total_groups_table <- cbind(group_1_table, group_2_table, group_3_table, group_4_table,
                            group_5_table, group_6_table, group_7_table, group_8_table)
head(total_groups_table)
colnames(total_groups_table) <- c("group_1", "group_2", "group_3", "group_4", "group_5",
                                  "group_6", "group_7", "group_8")
rownames(total_groups_table) <- c("gpos25", "gpos50", "gpos75", "gpos100")

group_1_table_perc <- group_1_table/sum(group_1_table)*100
group_2_table_perc <- group_2_table/sum(group_2_table)*100
group_3_table_perc <- group_3_table/sum(group_3_table)*100
group_4_table_perc <- group_4_table/sum(group_4_table)*100
group_5_table_perc <- group_5_table/sum(group_5_table)*100
group_6_table_perc <- group_6_table/sum(group_6_table)*100
group_7_table_perc <- group_7_table/sum(group_7_table)*100
group_8_table_perc <- group_8_table/sum(group_8_table)*100

total_groups_table_perc <- cbind(group_1_table_perc, group_2_table_perc, group_3_table_perc,
                                 group_4_table_perc,
                                 group_5_table_perc, group_6_table_perc, group_7_table_perc,
                                 group_8_table_perc)
head(total_groups_table_perc)
colnames(total_groups_table_perc) <- c("group_1", "group_2", "group_3", "group_4", "group_5",
                                       "group_6", "group_7", "group_8")
rownames(total_groups_table_perc) <- c("gpos25", "gpos50", "gpos75", "gpos100")

total_groups_table_perc <- round(total_groups_table_perc, 1)
