library(gtools)
library(maftools)
library(dplyr)
library(rstatix)

tcga.cohort = system.file('extdata', 'tcga_cohort.txt.gz', package = 'maftools')
tcga.cohort = data.table::fread(cmd = paste('zcat <', tcga.cohort), sep = '\t', stringsAsFactors = FALSE)
tcga.cohort$total <- tcga.cohort$total/50

source("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Code/Fig1/OS_by_quartile.R")
colnames(tcga.cohort)[2] <- "bcr_patient_barcode"
combinedDF <- merge(x = tcga.cohort, y = survival_data, by = 'bcr_patient_barcode')

combinedDF <- combinedDF[combinedDF$cohort %in% c("BRCA","UCEC","COAD","LGG","OV","SARC","THCA"), ]
combinedDF <- combinedDF[combinedDF$Quartiles %in% c(1,4), ]
combinedDF$Quartiles <- ifelse(combinedDF$Quartiles == 1, "Young", "Old")
combinedDF$Quartiles <- factor(combinedDF$Quartiles, levels = c("Young", "Old"))
combinedDF$cohort <- factor(combinedDF$cohort, levels = rev(c("THCA", "LGG", "BRCA", "SARC", "OV", "UCEC", "COAD")))


stat.test <- combinedDF %>%
  group_by(cohort) %>%
  wilcox_test(total ~ Quartiles) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(y.position = 3)
stat.test


p <- ggboxplot(data = combinedDF, x = 'Quartiles', y = 'total', fill = 'Quartiles', palette = "jco") +
  facet_wrap(~ cohort, nrow = 1) +
  scale_y_log10(name = "TMB", breaks = c(0.1,1,10,100,1000,10000),
                labels = c(1,2,3,4,5,6), limits = c(0.01,10000)) +
  labs(fill = 'Age') +
  theme_pubr(base_size =  35) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
        # axis.title.y = element_text(size = 20),
        # axis.text.y = element_text(size = 15),
        # strip.text = element_text(size = 20),
        # legend.title = element_text(size = 20),
        # legend.text =element_text(size = 15)
        )


stat.test$rounded <-  round(stat.test$p.adj, 10)
stat.test$fdrstar <- stars.pval(p.value = stat.test$rounded)
p + stat_pvalue_manual(stat.test, label = "fdrstar", hide.ns = T, label.size = 20)
# ggsave("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Figures/cancertype_tmb.jpg", dpi = 320, width = 12)
