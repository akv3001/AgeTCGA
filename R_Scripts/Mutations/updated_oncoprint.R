library(TCGAbiolinks)
library(maftools)
library(dplyr)

# Load data
cancer_code <- c("BRCA","COAD","THCA","UCEC","SARC","LGG","OV")
# Top 20 FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
          "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
          "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")
survival_data <- read.table("/Users/Yajas/Downloads/survival_data.txt", sep = "\t", row.names = 1,
                            stringsAsFactors = FALSE)
survival_data$Quartiles <- ifelse(survival_data$Quartiles == 1, "Young",
                                  ifelse(survival_data$Quartiles == 4, "Old", NA))
survival_data <- survival_data[!is.na(survival_data$Quartiles), ]
survival_data <- survival_data[survival_data$type %in% cancer_code, ]
colnames(survival_data)[2] <- "Tumor_Sample_Barcode"

# download MAF and combine
for (i in 1:length(cancer_code)) {
  if (i ==1){
    tumor_maf <- GDCquery_Maf(tumor = cancer_code[i],
                              pipelines = "mutect2")
  }
  if (i != 1){
    downloadedMAF <- GDCquery_Maf(tumor = cancer_code[i],
                                  pipelines = "mutect2")
    tumor_maf <- plyr::rbind.fill(tumor_maf, downloadedMAF)
    rm(downloadedMAF)
  }
}

# Create MAF file for samples of interest after removing FLAGS
tumor_maf$Tumor_Sample_Barcode <- substr(tumor_maf$Tumor_Sample_Barcode,1,12)
tumor_maf <- tumor_maf %>% 
  filter(!Hugo_Symbol %in% flags) %>% 
  filter(Tumor_Sample_Barcode %in% survival_data$Tumor_Sample_Barcode)
colnames(survival_data)[3] <- "Tumor Type"
colnames(survival_data)[35] <- "Age"
MAF <- read.maf(maf = tumor_maf, clinicalData = survival_data)

age_col <- setNames(c(pal_jco()(2)[1], pal_jco()(2)[2]), c("Young", "Old"))
cancer_col <- setNames(RColorBrewer::brewer.pal(n = 7,name = 'Dark2'), cancer_code)


# plot
par(mar=c(1,1,1,1))
# tiff("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Figures/updated_mutation/oncoplot.tiff",
#      width = 1000, height = 1000)
tiff("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Figures/updated_mutation/oncoplot.tiff",
     res = 300, width = 7, height = 7, units = "in")
oncoplot(maf =MAF, 
         top = 10, 
         bgCol = "white",
         legendFontSize = 1,sepwd_samples = 0,draw_titv = T,logColBar = T,
         sortByAnnotation = T,annotationFontSize = 1, fontSize = 1, titleFontSize = 1,
         annotationColor = list(Age = age_col, `Tumor Type` = cancer_col),
         gene_mar = 7, showTitle = FALSE,
         clinicalFeatures = c("Age","Tumor Type"), # which columns should be plotted
         annotationDat = survival_data) # data frame with metadata
dev.off()
