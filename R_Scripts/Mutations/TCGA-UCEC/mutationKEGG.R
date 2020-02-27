library(scales)
library(ggpubr)
library(viridis)
library(enrichplot)
library(DOSE)
library(ComplexHeatmap)
library(clusterProfiler)

data <- read.table("/Users/Yajas/Downloads/all_mutations_fdr_0.05.txt")
data <- data[data$cancercode == "UCEC", ]
data$Hugo_Symbol <- gsub('\\-.*','',data$Hugo_Symbol)

mart <- bitr(data$Hugo_Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Hugo_Symbol", "ENTREZ")
data <- merge(data, mart, "Hugo_Symbol")
data <- data[data$or > 1, ]

result <- enrichKEGG(use_internal_data = F, gene = data$ENTREZ, organism = "human", keyType = 'ncbi-geneid')
dotplot(result) + scale_color_viridis() + theme_bw(base_size = 25) +
  scale_y_discrete(labels = wrap_format(25))
# ggsave("/Users/Yajas/Downloads/UCEC_mutations_KEGG.jpg", height = 10)
