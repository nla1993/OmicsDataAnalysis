# Pie chart of new groups

library(ggplot2)
library(gridExtra)
total_groups_table_perc$gpos <- rownames(total_groups_table_perc)
# Source: Frequency table
G1 <- data.frame(total_groups_table_perc$gpos, total_groups_table_perc$group_1)
colnames(G1) <- c("class", "freq")
pie_1 <- ggplot(G1, aes(x = "", y=freq, fill = factor(class))) + 
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(round(freq), "%")), position = position_stack(vjust = 0.5)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Group 1") + 
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

pie_1 <- pie_1 + coord_polar(theta = "y") + theme(legend.position="none")

G2 <- data.frame(total_groups_table_perc$gpos, total_groups_table_perc$group_2)
colnames(G2) <- c("class", "freq")
pie_2 <- ggplot(G2, aes(x = "", y=freq, fill = factor(class))) + 
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(round(freq), "%")), position = position_stack(vjust = 0.5)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Group 2") +
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

pie_2 <- pie_2 + coord_polar(theta = "y") + theme(legend.position="none")

G3 <- data.frame(total_groups_table_perc$gpos, total_groups_table_perc$group_3)
colnames(G3) <- c("class", "freq")
pie_3 <- ggplot(G3, aes(x = "", y=freq, fill = factor(class))) + 
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(round(freq), "%")), position = position_stack(vjust = 0.5)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Group 3") +
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

pie_3 <- pie_3 + coord_polar(theta = "y") + theme(legend.position="none")

G4 <- data.frame(total_groups_table_perc$gpos, total_groups_table_perc$group_4)
colnames(G4) <- c("class", "freq")
pie_4 <- ggplot(G4, aes(x = "", y=freq, fill = factor(class))) + 
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(round(freq), "%")), position = position_stack(vjust = 0.5)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Group 4") +
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

pie_4 <- pie_4 + coord_polar(theta = "y") + theme(legend.position="none")


G5 <- data.frame(total_groups_table_perc$gpos, total_groups_table_perc$group_5)
colnames(G5) <- c("class", "freq")
pie_5 <- ggplot(G5, aes(x = "", y=freq, fill = factor(class))) + 
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(round(freq), "%")), position = position_stack(vjust = 0.5)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Group 5") +
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

pie_5<- pie_5 + coord_polar(theta = "y") + theme(legend.position="none")

G6 <- data.frame(total_groups_table_perc$gpos, total_groups_table_perc$group_6)
colnames(G6) <- c("class", "freq")
pie_6 <- ggplot(G6, aes(x = "", y=freq, fill = factor(class))) + 
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(round(freq), "%")), position = position_stack(vjust = 0.5)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Group 6") +
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))
pie_6<- pie_6 + coord_polar(theta = "y") + theme(legend.position="none")

G7 <- data.frame(total_groups_table_perc$gpos, total_groups_table_perc$group_7)
colnames(G7) <- c("class", "freq")
pie_7 <- ggplot(G7, aes(x = "", y=freq, fill = factor(class))) + 
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(round(freq), "%")), position = position_stack(vjust = 0.5)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Group 7") +
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

pie_7<- pie_7 + coord_polar(theta = "y") + theme(legend.position="none")

G8 <- data.frame(total_groups_table_perc$gpos, total_groups_table_perc$group_8)
colnames(G8) <- c("class", "freq")
pie_8 <- ggplot(G8, aes(x = "", y=freq, fill = factor(class))) + 
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = paste0(round(freq), "%")), position = position_stack(vjust = 0.5)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "Group 8") +
  scale_fill_manual(values = c("mediumpurple1", "chartreuse3", "coral2", "cyan3")) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

pie_8 <- pie_8 + coord_polar(theta = "y") #+ theme(legend.position="none")

grid.arrange(pie_1, pie_2, pie_3, pie_4, pie_5, pie_6, pie_7, pie_8, nrow = 2)

