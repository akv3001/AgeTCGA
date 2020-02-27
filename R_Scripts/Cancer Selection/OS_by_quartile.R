library(ggsci)
library(data.table)
library(openxlsx)
library(survival)
library(survminer)

survival_data <- read.xlsx("/Users/Yajas/Documents/Elemento/AgeTCGA-master/DATA/Immune Landscape of Cancer/TCGA-CDR-SupplementalTableS1.xltx")
survival_data <- survival_data[!is.na(survival_data$age_at_initial_pathologic_diagnosis), ]
survival_data$type[survival_data$type == "READ"] <- "COAD"

# create cancer type based age quartiles
setDT(survival_data)[,Quartiles := cut(age_at_initial_pathologic_diagnosis, quantile(age_at_initial_pathologic_diagnosis, probs = 0:4/4),
                                       labels = FALSE, include.lowest = TRUE), by = type]

# remove NA PFI instances/time
survival_data <- survival_data[!is.na(survival_data$OS),]
survival_data <- survival_data[!is.na(survival_data$OS.time),]

# loop to calculate OS survival
cancer_code <- unique(survival_data$type)
results_df <- data.frame(matrix(nrow = length(cancer_code),
                                ncol = 7))
colnames(results_df) <- c("Cancer", "Surv_N", "Surv_events", "hazard_ratio",
                          "lower_ci", "upper_ci", "Surv_pval")

for (i in 1:length(cancer_code)) {
  cancerDF <- survival_data[survival_data$type == cancer_code[i], ]
  
  surv_result <- summary(coxph(Surv(time= OS.time ,event = OS ) ~ Quartiles, 
                                  data = cancerDF))
  
  results_df$Surv_pval[i] <- coef(surv_result)[1,5]
  results_df$hazard_ratio[i] = round(surv_result$coefficients[2],5)
  results_df$Surv_N[i] <- surv_result$n
  results_df$Surv_events[i] <- surv_result$nevent
  results_df$lower_ci[i] <- surv_result$conf.int[3]
  results_df$upper_ci[i] <- surv_result$conf.int[4]
  results_df$Cancer[i] <- cancer_code[i]
}

# FDR Correct
results_df$Surv_pval <- p.adjust(results_df$Surv_pval, "fdr")
sig_cancers <- results_df$Cancer[results_df$Surv_pval<0.05]

results_df <- mutate(results_df, sig = ifelse(results_df$Surv_pval<0.05, "Sig.", "Not Sig."))

results_df <- results_df %>% arrange(sig, hazard_ratio)
levs <- results_df$Cancer
results_df$Cancer <- factor(results_df$Cancer, levels = levs)

ggplot(data = results_df, aes(x = Cancer, y = hazard_ratio, ymin = lower_ci, ymax = upper_ci)) + 
  geom_pointrange(aes(col = sig), shape = 18) +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, col = sig),width=0.5) +
  geom_hline(yintercept=1, lty=2) +
  scale_y_log10() +
  scale_color_aaas()+
  labs(y = "Hazard Ratio", x = "") +
  theme_pubr(legend = "none") + 
  coord_flip()
