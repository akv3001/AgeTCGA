library(pheatmap)
library(gplots)
library(viridis)
library(dplyr)
library(network)
library(ggplot2)
library(ggpubr)

#### for xcell

xcell_data <- read.csv("DATA/Deconvolution and ssGSEA/xCell_TCGA_RSEM.txt", stringsAsFactors = F, sep = "\t", header = T)    ## xcell tcga data
colnames(xcell_data) <-  gsub('\\.', '-', colnames(xcell_data))
filtered_meta <- Age_OnyCT.meta[which(Age_OnyCT.meta$CASE_ID %in% colnames(xcell_data)),]
length(which(colnames(xcell_data) %in% filtered_meta$CASE_ID))

##change rownames of xcell_data to the this:
##aDC Adipocytes Astrocytes B-cells Basophils CD4+ memory T-cells CD4+ naive T-cells 
##CD4+ T-cells CD4+ Tcm CD4+ Tem CD8+ naive T-cells CD8+ T-cells CD8+ Tcm CD8+ Tem cDC Chondrocytes 
##Class-switched memory B-cells CLP CMP DC Endothelial cells Eosinophils Epithelial cells Erythrocytes 
##Fibroblasts GMP Hepatocytes HSC iDC Keratinocytes ly Endothelial cells Macrophages Macrophages M1 Macrophages M2
##Mast cells Megakaryocytes Melanocytes Memory B-cells MEP Mesangial cells Monocytes 
##MPP MSC mv Endothelial cells Myocytes naive B-cells Neurons Neutrophils NK cells NKT Osteoblast pDC Pericytes
##Plasma cells Platelets Preadipocytes pro B-cells Sebocytes Skeletal muscle Smooth muscle Tgd cells Th1 cells Th2 cells Tregs


## different cell types
cell_types <- c("aDC","Adipocytes","Astrocytes","B-cells","Basophils",
                "CD4+ memory T-cells","CD4+ naive T-cells","CD4+ T-cells",
                "CD4+ Tcm","CD4+ Tem","CD8+ naive T-cells","CD8+ T-cells",
                "CD8+ Tcm","CD8+ Tem","cDC","Chondrocytes","Class-switched memory B-cells",
                "CLP","CMP","DC","Endothelial cells","Eosinophils","Epithelial cells",
                "Erythrocytes","Fibroblasts","GMP","Hepatocytes","HSC","iDC","Keratinocytes",
                "ly Endothelial cells","Macrophages","Macrophages M1","Macrophages M2",
                "Mast cells","Megakaryocytes","Melanocytes","Memory B-cells","MEP","Mesangial cells",
                "Monocytes","MPP","MSC","mv Endothelial cells","Myocytes","naive B-cells",
                "Neurons","Neutrophils","NK cells","NKT","Osteoblast","pDC","Pericytes",
                "Plasma cells","Platelets","Preadipocytes","pro B-cells","Sebocytes",
                "Skeletal muscle cells","Smooth muscle cells","Tgd cells","Th1 cells","Th2 cells","Tregs")

colnames(filtered_xcell) <- cell_types

## some more information to use if neededd
xcell_metadata <- as.data.frame(cbind(c("aDC","Adipocytes","Astrocytes","B-cells","Basophils","CD4+ memory T-cells","CD4+ naive T-cells","CD4+ T-cells","CD4+ Tcm","CD4+ Tem","CD8+ naive T-cells","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","cDC","Chondrocytes","Class-switched memory B-cells","CLP","CMP","DC","Endothelial cells","Eosinophils","Epithelial cells","Erythrocytes","Fibroblasts","GMP","Hepatocytes","HSC","iDC","Keratinocytes","ly Endothelial cells","Macrophages","Macrophages M1","Macrophages M2","Mast cells","Megakaryocytes","Melanocytes","Memory B-cells","MEP","Mesangial cells","Monocytes","MPP","MSC","mv Endothelial cells","Myocytes","naive B-cells","Neurons","Neutrophils","NK cells","NKT","Osteoblast","pDC","Pericytes","Plasma cells","Platelets","Preadipocytes","pro B-cells","Sebocytes","Skeletal muscle cells","Smooth muscle cells","Tgd cells","Th1 cells","Th2 cells","Tregs"),
                                c("Non-lymphocytes","Non-Hematopoietic","Non-Hematopoietic","Lymphocytes","Non-lymphocytes","Lymphocytes","Lymphocytes","Lymphocytes","Lymphocytes","Lymphocytes","Lymphocytes","Lymphocytes","Lymphocytes","Lymphocytes","Non-lymphocytes","Non-Hematopoietic","Lymphocytes","HSC","HSC","Non-lymphocytes","Non-Hematopoietic","Non-lymphocytes","Non-Hematopoietic","HSC","Non-Hematopoietic","HSC","Non-Hematopoietic","HSC","Non-lymphocytes","Non-Hematopoietic","Non-Hematopoietic","Non-lymphocytes","Non-lymphocytes","Non-lymphocytes","Non-lymphocytes","HSC","Non-Hematopoietic","Lymphocytes","HSC","Non-Hematopoietic","Non-lymphocytes","HSC","Non-Hematopoietic","Non-Hematopoietic","Non-Hematopoietic","Lymphocytes","Non-Hematopoietic","Non-lymphocytes","Lymphocytes","Lymphocytes","Non-Hematopoietic","Non-lymphocytes","Non-Hematopoietic","Lymphocytes","Non-lymphocytes","Non-Hematopoietic","Lymphocytes","Non-Hematopoietic","Non-Hematopoietic","Non-Hematopoietic","Lymphocytes","Lymphocytes","Lymphocytes","Lymphocytes"),
                                c("Myeloid","Stroma","Epithelial","Lymphoid","Myeloid","Lymphoid","Lymphoid","Lymphoid","Lymphoid","Lymphoid","Lymphoid","Lymphoid","Lymphoid","Lymphoid","Myeloid","Stroma","Lymphoid","HSC","HSC","Myeloid","Stroma","Myeloid","Epithelial","HSC","Stroma","HSC","Epithelial","HSC","Myeloid","Epithelial","Stroma","Myeloid","Myeloid","Myeloid","Myeloid","HSC","Epithelial","Lymphoid","HSC","Stroma","Myeloid","HSC","Stroma","Stroma","Stroma","Lymphoid","Epithelial","Myeloid","Lymphoid","Lymphoid","Stroma","Myeloid","Stroma","Lymphoid","HSC","Stroma","Lymphoid","Epithelial","Stroma","Stroma","Lymphoid","Lymphoid","Lymphoid","Lymphoid")), stringsAsFactors = F)


#### Look at subtype quartiles
levels(filtered_meta$TCGA_CancerCode)

## create data frames to store p values
results_df.p <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.p[,1] <- levels(filtered_meta$TCGA_CancerCode)
results_df.fdr <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.fdr[,1] <- levels(filtered_meta$TCGA_CancerCode)
results_df.mean <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.mean[,1] <- levels(filtered_meta$TCGA_CancerCode)

rownames(results_df) <- levels(filtered_meta$TCGA_CancerCode)

for (i in 1:length(levels(filtered_meta$TCGA_CancerCode))) {
  
  ## select case ids present in both datasets
  ## sort by case ids so everything matches up
  new_filtered_meta <- filtered_meta[filtered_meta$TCGA_CancerCode == levels(filtered_meta$TCGA_CancerCode)[i],]
  filtered_xcell <- xcell_data[,which(colnames(xcell_data) %in% new_filtered_meta$CASE_ID)]
  rownames(filtered_xcell) <- cell_types
  filtered_xcell <- filtered_xcell[,sort(colnames(filtered_xcell))]
  new_filtered_meta <- new_filtered_meta[order(new_filtered_meta$CASE_ID),]
  
  ## create a unified data frame
  compare_df <- data.frame(new_filtered_meta,t(filtered_xcell))
  
  ## setup for loop
  pval_list <- c()
  Q1mean.list <- c()
  Q3mean.list <- c()
  meanComp.list <- c()
  
  ## loop through each cell type and get p values
  ## loop from 14, look at colnames(compare_df)
  for (j in 14:ncol(compare_df)) {
    pval_list[j] <- wilcox.test(compare_df[which(compare_df$Subtype_Quartile == "Young_Q1"),j],
                                compare_df[which(compare_df$Subtype_Quartile == "Old_Q3"),j])$p.value
    Q1mean.list[j] <- mean(compare_df[which(compare_df$Subtype_Quartile == "Young_Q1"),j])
    Q3mean.list[j] <- mean(compare_df[which(compare_df$Subtype_Quartile == "Old_Q3"),j])
    meanComp.list[j] <- ifelse(Q1mean.list[j] > Q3mean.list[j], "< Quartile 1", 
                               ifelse(Q3mean.list[j] > Q1mean.list[j], "> Quartile 3", "Q1 =Q 3"))
  }
  
  ## create a data frame for p, fdr, etc
  pval_list <- pval_list[-c(1:13)]
  fdr_list <- p.adjust(pval_list, method = "fdr")
  mean_list <- meanComp.list[-c(1:13)]
  results_df.p[i,2:65] <- pval_list
  results_df.fdr[i,2:65] <- fdr_list
  results_df.mean[i,2:65] <- mean_list
}


colnames(results_df.p) <- c("TCGA_CancerCode", rownames(filtered_xcell))
colnames(results_df.fdr) <- c("TCGA_CancerCode", rownames(filtered_xcell))
colnames(results_df.mean) <- c("TCGA_CancerCode", rownames(filtered_xcell))
rownames(results_df.p) <- results_df.p[,1]
rownames(results_df.fdr) <- results_df.fdr[,1]
rownames(results_df.mean) <- results_df.mean[,1]

## heatmaps to visualize
pheatmap(as.matrix(t(results_df.p[2:ncol(results_df.p)])), main = "p")
pheatmap(as.matrix(t(results_df.fdr[2:ncol(results_df.fdr)])), main = "fdr")


## create a unified data frame for ggplot
dat0 =cbind(melt(results_df.mean,'TCGA_CancerCode'), melt(results_df.p,'TCGA_CancerCode')[,3], melt(results_df.fdr,'TCGA_CancerCode')[,3])
colnames(dat0) <- c("TCGA_CancerCode", "variable", "value", "p.val", "fdr")
dat0 <- dat0[which(dat0$fdr<0.05),]    ## only selecting significant values
ggplot(data = dat0, aes(x = TCGA_CancerCode, y = variable)) +
  geom_point(aes(size = p.val, shape = value, color = fdr)) +
  ggtitle("Subtype Quartile-based Age Stratification")+
  xlab("TCGA Cancer Code") +
  ylab("Cell Type")
ggsave("Figures/xCell TCGA Subtype Quartiles.jpg", dpi=320)
dev.off()

###### tertiles - no comments from here on but everything remains more or less the same
levels(filtered_meta$TCGA_CancerCode)
results_df.p <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.p[,1] <- levels(filtered_meta$TCGA_CancerCode)
results_df.fdr <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.fdr[,1] <- levels(filtered_meta$TCGA_CancerCode)
results_df.mean <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.mean[,1] <- levels(filtered_meta$TCGA_CancerCode)

rownames(results_df) <- levels(filtered_meta$TCGA_CancerCode)
for (i in 1:length(levels(filtered_meta$TCGA_CancerCode))) {
  new_filtered_meta <- filtered_meta[filtered_meta$TCGA_CancerCode == levels(filtered_meta$TCGA_CancerCode)[i],]
  filtered_xcell <- xcell_data[,which(colnames(xcell_data) %in% new_filtered_meta$CASE_ID)]
  rownames(filtered_xcell) <- cell_types
  filtered_xcell <- filtered_xcell[,sort(colnames(filtered_xcell))]
  new_filtered_meta <- new_filtered_meta[order(new_filtered_meta$CASE_ID),]
  compare_df <- data.frame(new_filtered_meta,t(filtered_xcell))
  
  pval_list <- c()
  Q1mean.list <- c()
  Q3mean.list <- c()
  meanComp.list <- c()
  for (j in 14:ncol(compare_df)) {
    pval_list[j] <- wilcox.test(compare_df[which(compare_df$Subtype_Tertile == "Young_T1"),j],
                                compare_df[which(compare_df$Subtype_Tertile == "Old_T3"),j])$p.value
    Q1mean.list[j] <- mean(compare_df[which(compare_df$Subtype_Tertile == "Young_T1"),j])
    Q3mean.list[j] <- mean(compare_df[which(compare_df$Subtype_Tertile == "Old_T3"),j])
    meanComp.list[j] <- ifelse(Q1mean.list[j] > Q3mean.list[j], "< T1", 
                               ifelse(Q3mean.list[j] > Q1mean.list[j], "> T2", "Tertile 1 = Tertile 3"))
  }
  pval_list <- pval_list[-c(1:13)]
  fdr_list <- p.adjust(pval_list, method = "fdr")
  mean_list <- meanComp.list[-c(1:13)]
  results_df.p[i,2:65] <- pval_list
  results_df.fdr[i,2:65] <- fdr_list
  results_df.mean[i,2:65] <- mean_list
}
colnames(results_df.p) <- c("TCGA_CancerCode", rownames(filtered_xcell))
colnames(results_df.fdr) <- c("TCGA_CancerCode", rownames(filtered_xcell))
colnames(results_df.mean) <- c("TCGA_CancerCode", rownames(filtered_xcell))
rownames(results_df.p) <- results_df.p[,1]
rownames(results_df.fdr) <- results_df.fdr[,1]
rownames(results_df.mean) <- results_df.mean[,1]

pheatmap(t(results_df.p[2:ncol(results_df.p)]), main = "p")
pheatmap(t(results_df.fdr[2:ncol(results_df.fdr)]), main = "fdr")


dat0 =cbind(melt(results_df.mean,'TCGA_CancerCode'), melt(results_df.p,'TCGA_CancerCode')[,3], melt(results_df.fdr,'TCGA_CancerCode')[,3])
colnames(dat0) <- c("TCGA_CancerCode", "variable", "value", "p.val", "fdr")
dat0 <- dat0[which(dat0$fdr<0.05),]
ggplot(data = dat0, aes(x = TCGA_CancerCode, y = variable)) + 
  geom_point(aes(size = p.val, shape = value, color = fdr)) +
  ggtitle("Subtype Tertile-based Age Stratification")+
  xlab("TCGA Cancer Code") +
  ylab("Cell Type")
ggsave("Figures/xCell TCGA Subtype Tertiles.jpg", dpi=320)
dev.off()


#### PanCan Quartiles
results_df.p <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.p[,1] <- levels(filtered_meta$TCGA_CancerCode)
results_df.fdr <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.fdr[,1] <- levels(filtered_meta$TCGA_CancerCode)
results_df.mean <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.mean[,1] <- levels(filtered_meta$TCGA_CancerCode)

filtered_xcell <- xcell_data[,which(colnames(xcell_data) %in% filtered_meta$CASE_ID)]
filtered_xcell <- filtered_xcell[,sort(colnames(filtered_xcell))]
filtered_meta <- filtered_meta[order(filtered_meta$CASE_ID),]
compare_df <- data.frame(filtered_meta,t(filtered_xcell))
colnames(compare_df)[14:ncol(compare_df)] <- cell_types

pval_list <- c()
Q1mean.list <- c()
Q3mean.list <- c()
meanComp.list <- c()

for (j in 14:ncol(compare_df)) {
  pval_list[j] <- wilcox.test(compare_df[which(compare_df$PanCan_Quartile == "Young_Q1"),j],
                              compare_df[which(compare_df$PanCan_Quartile == "Old_Q3"),j])$p.value
  Q1mean.list[j] <- mean(compare_df[which(compare_df$PanCan_Quartile == "Young_Q1"),j])
  Q3mean.list[j] <- mean(compare_df[which(compare_df$PanCan_Quartile == "Old_Q3"),j])
  meanComp.list[j] <- ifelse(Q1mean.list[j] > Q3mean.list[j], "< Q1", 
                             ifelse(Q3mean.list[j] > Q1mean.list[j], "> Q3", "Q1 = Q3"))
}
pval_list <- pval_list[-c(1:13)]
fdr_list <- p.adjust(pval_list, method = "fdr")
mean_list <- meanComp.list[-c(1:13)]

cc <- as.data.frame(cbind(colnames(compare_df)[14:ncol(compare_df)], mean_list, pval_list,fdr_list), stringsAsFactors = FALSE)
cc$pval_list <- as.numeric(cc$pval_list)
cc$fdr_list <- as.numeric(cc$fdr_list)

cc <- cc[which(cc$fdr_list < 0.05),]

ggplot(data = cc, aes(x = V1, y = -log10(pval_list))) + geom_point(aes(size = fdr_list, shape = mean_list))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Cell Type") + ylab("-lopg10 p value") + ggtitle("PanCan Quartile Stratification")
ggsave("Figures/xCell TCGA PanCan Quartiles.jpg", dpi=320)
dev.off()

#### PanCan Tertiles
results_df.p <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.p[,1] <- levels(filtered_meta$TCGA_CancerCode)
results_df.fdr <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.fdr[,1] <- levels(filtered_meta$TCGA_CancerCode)
results_df.mean <- data.frame(matrix(ncol = 65, nrow = length(levels(filtered_meta$TCGA_CancerCode))))
results_df.mean[,1] <- levels(filtered_meta$TCGA_CancerCode)

filtered_xcell <- xcell_data[,which(colnames(xcell_data) %in% filtered_meta$CASE_ID)]
filtered_xcell <- filtered_xcell[,sort(colnames(filtered_xcell))]
filtered_meta <- filtered_meta[order(filtered_meta$CASE_ID),]
compare_df <- data.frame(filtered_meta,t(filtered_xcell))
colnames(compare_df)[14:ncol(compare_df)] <- cell_types

pval_list <- c()
Q1mean.list <- c()
Q3mean.list <- c()
meanComp.list <- c()

for (j in 14:ncol(compare_df)) {
  pval_list[j] <- wilcox.test(compare_df[which(compare_df$PanCan_Tertile == "Young_T1"),j],
                              compare_df[which(compare_df$PanCan_Tertile == "Old_T3"),j])$p.value
  Q1mean.list[j] <- mean(compare_df[which(compare_df$PanCan_Tertile == "Young_T1"),j])
  Q3mean.list[j] <- mean(compare_df[which(compare_df$PanCan_Tertile == "Old_T3"),j])
  meanComp.list[j] <- ifelse(Q1mean.list[j] > Q3mean.list[j], "< T1", 
                             ifelse(Q3mean.list[j] > Q1mean.list[j], "> T2", "T1 = T3"))
}
pval_list <- pval_list[-c(1:13)]
fdr_list <- p.adjust(pval_list, method = "fdr")
mean_list <- meanComp.list[-c(1:13)]

cc <- as.data.frame(cbind(colnames(compare_df)[14:ncol(compare_df)], mean_list, pval_list,fdr_list), stringsAsFactors = FALSE)
cc$pval_list <- as.numeric(cc$pval_list)
cc$fdr_list <- as.numeric(cc$fdr_list)

cc <- cc[which(cc$fdr_list < 0.05),]

ggplot(data = cc, aes(x = V1, y = -log10(pval_list))) + geom_point(aes(size = fdr_list, shape = mean_list))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Cell Type") + ylab("-lopg10 p value") + ggtitle("PanCan Tertile Stratification")
ggsave("Figures/xCell TCGA PanCan Tertiles.jpg", dpi=320)
dev.off()
