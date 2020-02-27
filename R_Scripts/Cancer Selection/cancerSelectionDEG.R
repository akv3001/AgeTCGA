################################################################################################
# Load packages and data
################################################################################################

library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(limma)
library(edgeR)
library(TCGAutils)
library(biomaRt)
library(data.table)
library(org.Hs.eg.db)
library(ggrepel)
library(fgsea)
library(gridExtra)
library(fgsea)
library(pheatmap)
library(org.Hs.eg.db)
source("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Code/Fig1/OS_by_quartile.R")
load("~/Documents/Elemento/AgeTCGA-master/DATA/NewTCGAcounts/Full Matrix.RData")
pathways.hallmark <- gmtPathways("~/Documents/Elemento/AgeTCGA-master/DATA/GSEA_pathways/h.all.v7.0.symbols.gmt")

################################################################################################
# Unify datasets and select samples
################################################################################################

filtered_counts <- RNAseq_counts

# survival_data <- read.xlsx("/Users/Yajas/Documents/Elemento/AgeTCGA-master/DATA/Immune Landscape of Cancer/TCGA-CDR-SupplementalTableS1.xltx")
survival_data <- survival_data[!is.na(survival_data$age_at_initial_pathologic_diagnosis), ]
survival_data$type[survival_data$type == "READ"] <- "COAD"

# create cancer type based age quartiles
setDT(survival_data)[,Subtype_Quartile := cut(age_at_initial_pathologic_diagnosis, quantile(age_at_initial_pathologic_diagnosis, probs = 0:4/4),
                                              labels = FALSE, include.lowest = TRUE), by = type]
survival_data$Subtype_Quartile <- paste0("Q",survival_data$Subtype_Quartile)

# only look at primary solid tumors
filtered_counts <- filtered_counts[ ,(TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Solid Tumor" |
                                        TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Blood Derived Cancer - Peripheral Blood")]

# Remove duplicates if they exist
biospec_table <- TCGAbiospec(colnames(filtered_counts))
biospec_table <- biospec_table[order(biospec_table$plate,decreasing = TRUE), ]
filtered_counts <- filtered_counts[,colnames(filtered_counts)[order(biospec_table$plate, decreasing = TRUE)]]

if (sum(duplicated(biospec_table$submitter_id)) > 0){
  filtered_counts <- filtered_counts[, !duplicated(biospec_table$submitter_id)]
}

# change column names to match metadata
colnames(filtered_counts) <- substr(colnames(filtered_counts), 1, 12)

## selecting case IDs common to both datasets and discarding the rest
common_cases <- intersect(colnames(filtered_counts), survival_data$bcr_patient_barcode)
filtered_age <- survival_data[survival_data$bcr_patient_barcode %in% common_cases,]
filtered_counts <- filtered_counts[,colnames(filtered_counts) %in% common_cases]

## sorting case IDs
filtered_age <- filtered_age[order(filtered_age$bcr_patient_barcode),]
filtered_counts <- filtered_counts[,order(colnames(filtered_counts))]

rm("RNAseq_counts")


################################################################################################
# Create results dataframes
################################################################################################

## unique 3 letter codes
cancer_types <- levels(factor(filtered_age$type))
cancer_types <- sig_cancers
# cancer_types <- c("blca","brca","coad","hnsc","kirc","skcm","thca","ucec")

# This dataframe has the DE results
resultsDF_RNA <- as.data.frame(matrix(nrow = length(cancer_types), ncol = 10))
colnames(resultsDF_RNA) <- c("Q1_n","Q4_n","no_of_total_genes_tested", "no_of_sig_genes","sig_gene_fraction","Q1_min_age","Q1_max_age","Q4_min_age","Q4_max_age", "Cancer_Code")
rownames(resultsDF_RNA) <- cancer_types

# This dataframe has the FGSEA results
NES_DF <- as.data.frame(matrix(nrow = length(pathways.hallmark),
                               ncol = length(cancer_types)), 
                        row.names = names(pathways.hallmark))
colnames(NES_DF) <- cancer_types
NES_FDR_DF <- NES_DF


################################################################################################
# Differential Expression and GSEA
################################################################################################

## unique 3 letter codes
cancer_types <- levels(factor(filtered_age$type))
cancer_types <- sig_cancers
# cancer_types <- c("blca","brca","coad","hnsc","kirc","skcm","thca","ucec")
volcano_list <- list()
upregulated_genes <- list()
downregulated_genes <- list()

transcription_tt <- setNames(replicate(length(cancer_types),data.frame()),cancer_types)

for (i in 1:length(cancer_types)) {
  
  #### create a subset of the age and rna seq datasets
  #### select the current cancer type subjects which have which are either <Q1 or>Q3 wrt age
  age_subset <- filtered_age[which(filtered_age$type == cancer_types[i] & filtered_age$Subtype_Quartile != "Q2" &
                                     filtered_age$Subtype_Quartile != "Q3"), 
                             c("bcr_patient_barcode","Subtype_Quartile","age_at_initial_pathologic_diagnosis")]
  rna_subset <- filtered_counts[,which(colnames(filtered_counts) %in% age_subset$bcr_patient_barcode)]
  
  
  ### get min and max ages for each quartile
  resultsDF_RNA$Q1_min_age[i] <- min(as.numeric(age_subset$age_at_initial_pathologic_diagnosis[which(age_subset$Subtype_Quartile == "Q1")]))
  resultsDF_RNA$Q1_max_age[i] <- max(as.numeric(age_subset$age_at_initial_pathologic_diagnosis[which(age_subset$Subtype_Quartile == "Q1")]))
  resultsDF_RNA$Q4_min_age[i] <- min(as.numeric(age_subset$age_at_initial_pathologic_diagnosis[which(age_subset$Subtype_Quartile == "Q4")]))
  resultsDF_RNA$Q4_max_age[i] <- max(as.numeric(age_subset$age_at_initial_pathologic_diagnosis[which(age_subset$Subtype_Quartile == "Q4")]))
  
  
  #### limma
  input <- as.data.frame(age_subset[,c("bcr_patient_barcode","Subtype_Quartile")])
  input$Subtype_Quartile <- factor(input$Subtype_Quartile, levels = c("Q4","Q1"))
  condition <- c("Q4","Q1")
  comparison <- "Q4-Q1"
  design <- model.matrix(~0+ input$Subtype_Quartile)
  colnames(design) <- condition
  contmatrix <- makeContrasts(as.character(comparison),levels=design)
  resultsDF_RNA$Q4_n[i] <- table(design)[1]
  resultsDF_RNA$Q1_n[i] <- table(design)[2]
  
  Count_input <- rna_subset
  Count_input<-Count_input[which(rowSums(cpm(Count_input)> 2) >= 2) ,]
  resultsDF_RNA$no_of_total_genes_tested[i] <- nrow(rna_subset)
  print(dim(Count_input))
  Data.voom <- voom(Count_input,plot=FALSE)
  Voom_output <- Data.voom
  fit <- lmFit(Data.voom,design)
  fit2 <- contrasts.fit(fit,contmatrix)
  fit2 <-eBayes(fit2)
  tt <- topTable(fit2,n=Inf)
  resultsDF_RNA$no_of_sig_genes[i] <- length(which(tt$adj.P.Val<0.05))
  resultsDF_RNA$sig_gene_fraction[i] <- resultsDF_RNA$no_of_sig_genes[i]/resultsDF_RNA$no_of_total_genes_tested[i]
  
  # Run fgsea
  set.seed(2019)
  tt[,7] <- mapIds(org.Hs.eg.db, keys = rownames(tt), keytype = "ENSEMBL", column="SYMBOL", "first")
  tt <- na.omit(tt)
  tt <- tt[-which(duplicated(tt$V7)),]
  rownames(tt) <- tt$V7
  
  # Pairwise analysis
  transcription_tt[[i]] <- tt
  
  tt <- tt[order(tt$t, decreasing = T),]
  ranks_TCGA <- tt$t
  names(ranks_TCGA) <- rownames(tt)
  fgseaRes_TCGA <- fgsea(pathways=pathways.hallmark, stats=ranks_TCGA, nperm=1000)
  NES_FDR_DF[,i] <- fgseaRes_TCGA$padj
  NES_DF[,i] <- ifelse(fgseaRes_TCGA$padj<0.05, fgseaRes_TCGA$NES,0)
  
  resultsDF_RNA$Cancer_Code[i] <- cancer_types[i]
  
  
  upregulated_genes[[i]] <- ifelse(tt$logFC>0 & tt$adj.P.Val<0.05, rownames(tt)[tt$logFC>0 & tt$adj.P.Val<0.05], NA)
  downregulated_genes[[i]] <- ifelse(tt$logFC<0 & tt$adj.P.Val<0.05, rownames(tt)[tt$logFC<0 & tt$adj.P.Val<0.05], NA)
  
  # volcano plot
  tt <- mutate(tt, neglogfdr = -log10(adj.P.Val))
  
  tt <- mutate(tt, cutoffs = ifelse(adj.P.Val < 0.05 & abs(logFC) > 0.5, "p & FDR",
                                    ifelse(adj.P.Val < 0.05 & abs(logFC) <= 0.5, "p",
                                           ifelse(adj.P.Val >= 0.05 & abs(logFC) > 0.5, "FDR", "none"))))
  tt <- mutate(tt,cancer = cancer_types[i])
  volcano_list[i] <- ggplot(tt, aes(x = logFC, y = neglogfdr)) + geom_point()
  
}

names(volcano_list) <- cancer_types
names(upregulated_genes) <- cancer_types
names(downregulated_genes) <- cancer_types


################################################################################################
# Inspect Results
################################################################################################

head(resultsDF_RNA)
rownames(NES_DF) <- substring(rownames(NES_DF),10)
head(NES_DF)


################################################################################################
# DE Barplot and subset cancer types of interest
################################################################################################

resultsDF_RNA$sig_gene_fraction <- resultsDF_RNA$sig_gene_fraction*100
resultsDF_RNA <- mutate(resultsDF_RNA, Phenotype = ifelse(resultsDF_RNA$sig_gene_fraction > 1,
                                                          "Age-associated","Not age-associated"))

plotOrder <- resultsDF_RNA$Cancer_Code[order(resultsDF_RNA$no_of_sig_genes)]

rownames(resultsDF_RNA) <- resultsDF_RNA$Cancer_Code
resultsDF_RNA <- resultsDF_RNA[plotOrder, ]

ggbarplot(resultsDF_RNA, x = 'Cancer_Code', y = 'sig_gene_fraction',
          x.text.angle = 90, fill = 'Phenotype', title = "Differential Expression",
          ylab = "Percent DEG", xlab = F, color = pal_aaas()(1)) + 
  geom_hline(yintercept = 1, lty = 2, lwd = 2) +
  scale_fill_manual(values = pal_jama()(4)[3:4]) +
  theme_pubr(base_size = 35) +
  theme(# axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    legend.direction = "vertical",
    legend.justification = c("left", "top"),
    legend.position = c(.3, .35)) +
  coord_flip()
# ggsave("Final/Figures/Figure 1/Percent_DEG.jpg", dpi = 320)
# ggsave("Final/Figures/Figure 1/Percent_DEG.jpg", dpi = 320, width = 9, height = 9)

plotOrder <- resultsDF_RNA$Cancer_Code[order(resultsDF_RNA$no_of_sig_genes)]

c.types <- resultsDF_RNA$Cancer_Code
c.types <- plotOrder[8:length(plotOrder)]
c.types <- rownames(resultsDF_RNA)[resultsDF_RNA$sig_gene_fraction>1]


################################################################################################
# GSEA Heatmap
################################################################################################

hallmark_annotations <- read.csv("/Users/Yajas/Documents/Elemento/AgeTCGA-master/DATA/GSEA_pathways/HALLMARK_ANNOTATIONS.csv", 
                                 header = T, stringsAsFactors = F, sep = ",")
NES_DF <- NES_DF[,colnames(NES_DF) %in% c.types]
NES_DF <- NES_DF[order(rownames(NES_DF)),]
hallmark_annotations <- hallmark_annotations[order(hallmark_annotations$HALLMARK_PATHWAY),]
all(rownames(NES_DF) == hallmark_annotations$HALLMARK_PATHWAY)

keepPathways <- c(which(hallmark_annotations$ANNOTATION == "IMMUNE"),which(hallmark_annotations$ANNOTATION == "PROLIFERATION"),which(hallmark_annotations$ANNOTATION == "DNA DAMAGE"))

NES_DF <- NES_DF[keepPathways, ]
hallmark_annotations <- hallmark_annotations[keepPathways,]
row_anno <- as.data.frame(matrix(hallmark_annotations$ANNOTATION),row.names = hallmark_annotations$HALLMARK_PATHWAY)
colnames(row_anno) <- "Pathways"

# Set annotation colors
library(RColorBrewer)
library(ggsci)
mycolors <- list(`Pathways` = c(`DEVELOPMENT` = brewer.pal(8, "Set3")[1],
                                `IMMUNE` = brewer.pal(8, "Set3")[2],
                                `SIGNALING` = brewer.pal(8, "Set3")[3],
                                `CELLULAR COMPONENT` = brewer.pal(8, "Set3")[4],
                                `PATHWAY` = brewer.pal(8, "Set3")[5],
                                `METABOLIC` = brewer.pal(8, "Set3")[6],
                                `DNA DAMAGE` = brewer.pal(8, "Set3")[7],
                                `PROLIFERATION` = brewer.pal(8, "Set3")[8]),
                 `High DEG` = c(Yes = "dodgerblue2", No = "white"))


rownames(NES_DF) <- gsub('\\_', ' ', rownames(NES_DF))
tiff("Final/Figures/Figure 1/Immune_heatmap.tiff", res = 320, height = 11, width = 15, units = "in")
pheatmap(mat = NES_DF, 
         border_color = "white", cellwidth = 35,cellheight = 35, fontsize = 35,
         color = pal_gsea("default")(12), 
         show_rownames = T, show_colnames = T,annotation_names_row = F, cutree_rows = 2, legend_breaks = c(-2,-1,0,1,2,3),
         legend_labels = c("Young", " ", "","", "", "Old"))
dev.off()


################################################################################################
# Volcano Plots
################################################################################################

df <- bind_rows(volcano_list[c.types])
df$neglogpval <- -log10(df$P.Value)
df$cancer <- factor(df$cancer, levels = plotOrder)
# df <- mutate(df, drivers = df$V7 %in% driverM)
df$cutoffs <- ifelse(df$logFC > 0 & df$adj.P.Val < 0.05, "Upregulated",
                     ifelse(df$logFC < 0 & df$adj.P.Val < 0.05, "Downregulated", "No change"))

volc <- ggscatter(
  df,
  x = 'logFC',
  y = 'neglogfdr',
  size = 0.2,
  alpha = 1/3,
  xlab = FALSE,
  ylab = expression(~-log[10]~FDR),
  color = 'cutoffs',
  palette = c("dodgerblue2", "black", "firebrick2")
) +
  theme_pubr(base_size = 35) +
  facet_wrap('cancer', nrow = 1, scales = "free_y") +
  theme(
    legend.position = "none",
    strip.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_x_continuous(breaks = c(-2,0,2), labels = c(-2,0,2))


metaDF <-
  survival_data[survival_data$type %in% c.types &
                  survival_data$Subtype_Quartile != "Q2" & 
                  survival_data$Subtype_Quartile != "Q3",]
metaDF$Subtype_Quartile <- ifelse(metaDF$Subtype_Quartile == "Q1", "Young", "Old")
metaDF$Subtype_Quartile <- factor(metaDF$Subtype_Quartile, levels = c("Young", "Old"))
metaDF$type <- factor(metaDF$type, levels = plotOrder)

boxplots <-
  ggboxplot(
    metaDF,
    x = 'Subtype_Quartile',
    y = 'age_at_initial_pathologic_diagnosis',
    xlab = FALSE,
    ylab = "Age",
    fill = 'Subtype_Quartile',
    palette = 'jco',
    ylim = c(10, 100),
  ) +
  theme_pubr(base_size = 35) +
  facet_wrap('type', nrow = 1) + 
  labs(fill = "Age Group") +
  theme(
    legend.position = "top",
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_y_continuous(breaks = c(25, 75), labels = c(25, 75))


ggarrange(boxplots, volc, nrow = 2,heights = c(1,2),align = "v") # ggsave("Final/Figures/Figure 1/Volcano.jpg", dpi = 320, width = 30, height = 11)
