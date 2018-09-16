# Abundance of H1 variants at new groups
###########################
head(group_1)
boxplot(group_1$H1.2_HeLa, group_1$H1X_HeLa,
        group_2$H1.2_HeLa, group_2$H1X_HeLa,
        group_3$H1.2_HeLa, group_3$H1X_HeLa,
        group_4$H1.2_HeLa, group_4$H1X_HeLa,
        group_5$H1.2_HeLa, group_5$H1X_HeLa,
        group_6$H1.2_HeLa, group_6$H1X_HeLa,
        group_7$H1.2_HeLa, group_7$H1X_HeLa,
        group_8$H1.2_HeLa, group_8$H1X_HeLa,
        col = c("sandybrown", "lightblue1","sandybrown", "lightblue1","sandybrown", "lightblue1",
                "sandybrown", "lightblue1","sandybrown", "lightblue1","sandybrown", "lightblue1",
                "sandybrown", "lightblue1","sandybrown", "lightblue1"),
        outline = F, main = "H1 variants abundance at new groups", xlab=NULL)
abline(h = 0, col = "grey", lwd = 1)
abline(h = 0.5, col = "forestgreen", lwd = 1)
abline(h = -0.5, col = "lightblue", lwd = 1)
abline(v=2.5, col = "grey", lwd = 1)
abline(v=4.5, col = "grey", lwd = 1)
abline(v=6.5, col = "grey", lwd = 1)
abline(v=8.5, col = "grey", lwd = 1)
abline(v=10.5, col = "grey", lwd = 1)
abline(v=12.5, col = "grey", lwd = 1)
abline(v=14.5, col = "grey", lwd = 1)


