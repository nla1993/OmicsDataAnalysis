# Scatter plots to see correlations

# First run scripts 1, 3 and 4
source("~/Documents/OneDrive - Universitat de Vic/TFM/Scripts_final/H1atCytobands/CytobandAggregation/1.ReadingFiles_PTMs&CRs_NORM.R")
source("~/Documents/OneDrive - Universitat de Vic/TFM/Scripts_final/H1atCytobands/CytobandAggregation/3.ReadingH1AbundanceFiles.R")
source("~/Documents/OneDrive - Universitat de Vic/TFM/Scripts_final/H1atCytobands/CytobandAggregation/4.PTMs_CRs_heatmap_NORM.R")

library(ggplot2)
library(gridExtra)

# merge PTMs_CRs and H1variants
H1variants_PTMs_CRs <- merge(PTMs_CRs_NORM, H1variants, by="row.names")
head(H1variants_PTMs_CRs)
rownames(H1variants_PTMs_CRs) <- H1variants_PTMs_CRs$Row.names
H1variants_PTMs_CRs <- H1variants_PTMs_CRs[, -1]
head(H1variants_PTMs_CRs)

# re-ordering
H1variants_PTMs_CRs$V1 <- rownames(H1variants_PTMs_CRs)
head(H1variants_PTMs_CRs)
order_3 <- rownames(PTMs_CRs_NORM)

H1variants_PTMs_CRs <- H1variants_PTMs_CRs[match(order_3, H1variants_PTMs_CRs$V1),]
head(H1variants_PTMs_CRs)
H1variants_PTMs_CRs <- H1variants_PTMs_CRs[, -15]

# scale the data
H1variants_PTMs_CRs_scaled <-  as.data.frame(apply(H1variants_PTMs_CRs, 2, scale))
head(H1variants_PTMs_CRs_scaled)
rownames(H1variants_PTMs_CRs_scaled) <- rownames(H1variants_PTMs_CRs)
head(H1variants_PTMs_CRs_scaled)

head(H1variants_PTMs_CRs_scaled)
scatter_hela <- H1variants_PTMs_CRs_scaled
scatter_hela$cytoband_group <- c(rep("gpos25", 86), rep("gpos50", 121), rep("gpos75", 89), rep("gpos100", 81))
head(scatter_hela)

#correlation tests

a <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$H3K9me3, method= "pearson")
b <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$H3K27me3, method= "pearson")
c <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$P300, method= "pearson")
d <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$RNAPOLII, method= "pearson")
i <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$H3K36me3, method= "pearson")
k <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$CTCF, method="pearson")
m <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$H3K4me3, method="pearson")
o <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$H3K27ac, method="pearson")
q <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$H3K9ac, method="pearson")
s <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$CHD1, method="pearson")
u <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$H3K4me1, method="pearson")
w <- cor.test(scatter_hela$H1.2_HeLa, scatter_hela$EZH2, method="pearson")


e <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$H3K9me3, method= "pearson")
f <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$H3K27me3, method= "pearson")
g <-cor.test(scatter_hela$H1X_HeLa, scatter_hela$P300, method= "pearson")
h <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$RNAPOLII, method= "pearson")
j <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$H3K36me3, method= "pearson")
l <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$CTCF, method="pearson")
n <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$H3K4me3, method="pearson")
p <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$H3K27ac, method="pearson")
r <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$H3K9ac, method="pearson")
t <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$CHD1, method="pearson")
v <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$H3K4me1, method="pearson")
x <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$EZH2, method="pearson")


y <- cor.test(scatter_hela$H1X_HeLa, scatter_hela$H1.2_HeLa, method= "pearson")

corr_table <- data.frame(x= c(a$p.value, b$p.value, c$p.value, d$p.value, i$p.value,
                              k$p.value, m$p.value, o$p.value, q$p.value, s$p.value,
                              u$p.value, w$p.value,
                              e$p.value, f$p.value, g$p.value, h$p.value, j$p.value,
                              l$p.value, n$p.value, p$p.value, r$p.value, t$p.value,
                              v$p.value, x$p.value,
                              y$p.value),
                         y= c(a$estimate, b$estimate, c$estimate, d$estimate, i$estimate,
                              k$estimate, m$estimate, o$estimate, q$estimate, s$estimate,
                              u$estimate, w$estimate,
                              e$estimate, f$estimate, g$estimate, h$estimate, j$estimate,
                              l$estimate, n$estimate, p$estimate, r$estimate, t$estimate,
                              v$estimate, x$estimate,
                              y$estimate),
                         z= c(rep("H1.2", 12), rep("H1X", 13)))
corr_table <- as.data.frame(t(corr_table))
row.names(corr_table) <- c("p.value", "corr.coef", "VS.Variant")
colnames(corr_table) <- c("H3K9me3", "H3K27me3", "P300", "RNAPOL_II", "H3K36me3",
                          "CTCF", "H3K4me3", "H3K27ac", "H3K9ac", "CHD1", "H3K4me1", "EZH2",
                          "H3K9me3", "H3K27me3", "P300", "RNAPOL_II", "H3K36me3", 
                          "CTCF", "H3K4me3", "H3K27ac", "H3K9ac", "CHD1", "H3K4me1", "EZH2",
                          "H1X")

# qplot(RNAPOLII, H1.2_HeLa, colour = cytoband_group, data = scatter_hela) +
#   geom_smooth(method = "lm", se = FALSE)


H3K9me3_H1.2 <- ggplot(scatter_hela, aes(x=H3K9me3, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K9me3 vs H1.2", subtitle = "p.value = 1.52e-10 \ncorr.coef= 0.32")

H3K9me3_H1.2 <-H3K9me3_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K9me3_H1X <- ggplot(scatter_hela, aes(x=H3K9me3, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K9me3 vs H1X", subtitle = "p.value = 4.14e-5 \ncorr.coef= 0.21")

H3K9me3_H1X <-H3K9me3_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K27me3_H1.2 <- ggplot(scatter_hela, aes(x=H3K27me3, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K27me3 vs H1.2", subtitle = "p.value = 1.37e-1 \ncorr.coef= -0.07")

H3K27me3_H1.2 <- H3K27me3_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K27me3_H1X <- ggplot(scatter_hela, aes(x=H3K27me3, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "H3K27me3 vs H1X", subtitle = "p.value = 4.75e-2 \ncorr.coef= 0.11")

H3K27me3_H1X <-H3K27me3_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

P300_H1.2 <- ggplot(scatter_hela, aes(x=P300, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "P300 vs H1.2", subtitle = "p.value = 7.5e-3 \ncorr.coef= -0.13")

P300_H1.2 <- P300_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

P300_H1X <- ggplot(scatter_hela, aes(x=P300, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "P300 vs H1X", subtitle = "p.value = 1.99e-8 \ncorr.coef= -0.28")

P300_H1X <-P300_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

RNAPOLII_H1.2 <- ggplot(scatter_hela, aes(x=RNAPOLII, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "RNA POL II vs H1.2", subtitle = "p.value = 4.1e-10 \ncorr.coef= -0.31")

RNAPOLII_H1.2 <- RNAPOLII_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

RNAPOLII_H1X <- ggplot(scatter_hela, aes(x=RNAPOLII, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "RNA POL II vs H1X", subtitle = "p.value = 3.29e-01 \ncorr.coef= 0.05")

RNAPOLII_H1X <- RNAPOLII_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K36me3_H1.2 <- ggplot(scatter_hela, aes(x=H3K36me3, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "H3K36me3 vs H1.2", subtitle = "p.value = 2.01e-4  \ncorr.coef= -0.19")

H3K36me3_H1.2 <-H3K36me3_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K36me3_H1X <- ggplot(scatter_hela, aes(x=H3K36me3, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K36me3 vs H1X", subtitle = "p.value = 7.17e-1  \ncorr.coef= -0.02")

H3K36me3_H1X <-H3K36me3_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

CTCF_H1.2 <- ggplot(scatter_hela, aes(x=CTCF, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "CTCF vs H1.2", subtitle = "p.value = 9.86e-2  \ncorr.coef= -0.085")

CTCF_H1.2 <-CTCF_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

CTCF_H1X <- ggplot(scatter_hela, aes(x=CTCF, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "CTCF vs H1X", subtitle = "p.value = 3.65e-1  \ncorr.coef= 0.0467")

CTCF_H1X <-CTCF_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K4me3_H1.2 <- ggplot(scatter_hela, aes(x=H3K4me3, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K4me3 vs H1.2", subtitle = "p.value = 1.51e-07  \ncorr.coef= -0.266")

H3K4me3_H1.2 <-H3K4me3_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K4me3_H1X <- ggplot(scatter_hela, aes(x=H3K4me3, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K4me3 vs H1X", subtitle = "p.value = 0.157  \ncorr.coef= -0.073")

H3K4me3_H1X <-H3K4me3_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K27ac_H1.2 <- ggplot(scatter_hela, aes(x=H3K27ac, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "H3K27ac vs H1.2", subtitle = "p.value = 7.34e-09  \ncorr.coef= -0.292")

H3K27ac_H1.2 <-H3K27ac_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K27ac_H1X <- ggplot(scatter_hela, aes(x=H3K27ac, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K27ac vs H1X", subtitle = "p.value = 0.0258  \ncorr.coef= -0.114")

H3K27ac_H1X <-H3K27ac_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K9ac_H1.2 <- ggplot(scatter_hela, aes(x=H3K9ac, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K9ac vs H1.2", subtitle = "p.value = 5.29e-08  \ncorr.coef= -0.2752")

H3K9ac_H1.2 <-H3K9ac_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))


H3K9ac_H1X <- ggplot(scatter_hela, aes(x=H3K9ac, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "H3K9ac vs H1X", subtitle = "p.value = 0.952  \ncorr.coef= -0.003")

H3K9ac_H1X <-H3K9ac_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

CHD1_H1.2 <- ggplot(scatter_hela, aes(x=CHD1, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "CHD1 vs H1.2", subtitle = "p.value = 3.13e-09  \ncorr.coef= -0.299")

CHD1_H1.2 <-CHD1_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

CHD1_H1X <- ggplot(scatter_hela, aes(x=CHD1, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "CHD1 vs H1X", subtitle = "p.value = 0.315  \ncorr.coef= 0.051")

CHD1_H1X <-CHD1_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K4me1_H1.2 <- ggplot(scatter_hela, aes(x=H3K4me1, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K4me1 vs H1.2", subtitle = "p.value = 6.68e-07  \ncorr.coef= -0.252")

H3K4me1_H1.2 <-H3K4me1_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

H3K4me1_H1X <- ggplot(scatter_hela, aes(x=H3K4me1, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H3K4me1 vs H1X", subtitle = "p.value = 0.002  \ncorr.coef= -0.154")

H3K4me1_H1X <- H3K4me1_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

EZH2_H1.2 <- ggplot(scatter_hela, aes(x=EZH2, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "EZH2 vs H1.2", subtitle = "p.value = 0.015  \ncorr.coef= 0.124")

EZH2_H1.2<- EZH2_H1.2 + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

EZH2_H1X <- ggplot(scatter_hela, aes(x=EZH2, y=H1X_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = F,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) + 
  ggtitle(label = "EZH2 vs H1X", subtitle = "p.value = 0.07  \ncorr.coef= 0.093")

EZH2_H1X <-EZH2_H1X + scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

#########################################################################
ggplot(scatter_hela, aes(x=H1X_HeLa, y=H1.2_HeLa)) + 
  geom_smooth(method = "lm", se = FALSE, color="black") + 
  geom_point(aes(colour=cytoband_group),  show.legend = T,              # colour depends on cond2
             size=1) + theme_classic() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(), plot.title = element_text(hjust = 0.5, color = "#666666"),
        plot.subtitle = element_text(hjust = 0.5, color = "#666666")) +
  ggtitle(label = "H1.2 vs H1X", subtitle = "p.value = 1.23e-62 \ncorr.coef= 0.72") + 
  scale_color_manual(values=c("mediumpurple1", "chartreuse3", "coral2", "cyan3"))

grid.arrange(H3K9me3_H1.2, H3K27me3_H1.2, H3K4me3_H1.2, H3K4me1_H1.2, H3K27ac_H1.2, H3K9ac_H1.2,
            nrow=2)
grid.arrange(H3K36me3_H1.2, EZH2_H1.2, CTCF_H1.2, P300_H1.2, CHD1_H1.2, RNAPOLII_H1.2, nrow=2)
grid.arrange(H3K9me3_H1X, H3K27me3_H1X, H3K4me3_H1X, H3K4me1_H1X, H3K27ac_H1X, H3K9ac_H1X, nrow=2)
grid.arrange(H3K36me3_H1X, EZH2_H1X, CTCF_H1X, P300_H1X, CHD1_H1X, RNAPOLII_H1X,  nrow=2)

