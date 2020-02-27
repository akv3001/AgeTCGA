library(ggsci)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(gridExtra)

young_colour <- "#0073C2FF"    # jco blue
old_colour <- "#EFC000FF"    # jco yellow

# normal curve to show young/old
normal_plot <- ggplot(NULL, aes(c(-2.5,2.5))) +
  geom_area(stat = "function", fun = dnorm, fill = young_colour,color = young_colour, xlim = c(-3, -1)) + 
  geom_area(stat = "function", fun = dnorm, fill = "white", xlim = c(-1, 1)) +
  geom_area(stat = "function", fun = dnorm, fill = old_colour,color = old_colour, xlim = c(1, 3)) +
  geom_line(stat = "function", fun = dnorm,,color = "black", xlim = c(-3, -1)) +
  geom_line(stat = "function", fun = dnorm, color = "black", xlim = c(-1, 1)) +
  geom_line(stat = "function", fun = dnorm, color = "black", xlim = c(1, 3)) +
  labs(x = "TCGA-Cancer Age", y = "") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = c(-2,2), labels = c("Young", "Old")) + 
  theme_pubr(base_size = 35) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(face = "italic"))

# Setup 2x2 table to show age
ageTable <- data.frame(matrix(nrow = 3, ncol = 2))
dimnames(ageTable) <- list(c("TCGA", "METABRIC","GTEx"), c("Young","Old"))
ageTable$Young <- c("Quartile 1", "TCGA-BRCA \n Quartile 1", "< 50 years")
ageTable$Old <- c("Quartile 4",  "TCGA-BRCA \n Quartile 4", "> 59 years")


## Quartile plot
# load("~/Documents/Elemento/AgeTCGA-master/DATA/cBioportal_Survival_Data/age_metadata.RData")
source("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Code/Fig1/OS_by_quartile.R")

cancer_code <- unique(survival_data$type)

quartiles_DF <- survival_data %>%
  group_by(type) %>%
  summarize(q1_min = min(age_at_initial_pathologic_diagnosis[Quartiles == 1]),
            q1_max = max(age_at_initial_pathologic_diagnosis[Quartiles == 1]),
            q3_min = min(age_at_initial_pathologic_diagnosis[Quartiles == 4]),
            q3_max = max(age_at_initial_pathologic_diagnosis[Quartiles == 4]),
            median = median(age_at_initial_pathologic_diagnosis))

# order the dataframe
quartiles_DF$type <- factor(quartiles_DF$type,
                                       levels = quartiles_DF$type[order(quartiles_DF$median)])
# plot quartiles
quartilePlot <- ggplot(quartiles_DF) + 
  geom_rect(aes(xmin = type,
                xmax = type,
                ymin = q1_min,
                ymax = q1_max),
            color = young_colour,
            size = 1) +
  geom_rect(aes(xmin = type,
                xmax = type,
                ymin = q3_min,
                ymax = q3_max),
            color = old_colour,
            size = 1)+
  geom_point(aes(x = type, y = median)) +
  coord_flip() +
  theme_pubr() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  labs(x = "TCGA \n Cancer Type", y = "Age")


Tdefault <- ttheme_default(base_size = 11.5)
normal <- normal_plot + annotation_custom(tableGrob(ageTable, theme = Tdefault),
                                xmin = -8, ymin = 0.25) 

library(cowplot)
ggdraw(normal) + draw_plot(quartilePlot, x = 0.32, y = 0.3, scale = 0.3)


ggsave("Final/Figures/Figure 1/Combined_Schematic.jpg", dpi = 320)
