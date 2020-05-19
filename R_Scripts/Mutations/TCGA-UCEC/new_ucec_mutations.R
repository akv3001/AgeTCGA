library(gtable)
library(grid)
library(gtools)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(TCGAbiolinks)
library(maftools)

drivers <- readxl::read_xlsx("/Users/Yajas/Documents/Elemento/AgeTCGA-master/DATA/Immune Landscape of Cancer/driver.genes.xlsx",
                             sheet = 2, trim_ws = TRUE, skip = 3)[,1:2]
survival_data <- read.table("/Users/Yajas/Downloads/survival_data.txt", sep = "\t", row.names = 1,
                            stringsAsFactors = FALSE)
survival_data <- survival_data %>%
  dplyr::filter(type == "UCEC") %>%
  dplyr::filter(Quartiles %in% c(1,4)) %>%
  dplyr::mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"))
colnames(survival_data)[2] <- "Tumor_Sample_Barcode"

tumor_maf <- GDCquery_Maf(tumor = "UCEC",
                          pipelines = "mutect2")
tumor_maf$Tumor_Sample_Barcode <- substr(tumor_maf$Tumor_Sample_Barcode,1,12)


# Get young and old barcodes
young_samples <- survival_data$Tumor_Sample_Barcode[survival_data$Quartiles == "Young"]
old_samples <- survival_data$Tumor_Sample_Barcode[survival_data$Quartiles == "Old"]

# clinical subsets
young_clinical <- survival_data %>% filter(Quartiles == "Young")
old_clinical <- survival_data %>% filter(Quartiles == "Old")

# maf subsets
young_maf <- tumor_maf %>% filter(Tumor_Sample_Barcode %in% young_samples)
old_maf <- tumor_maf %>% filter(Tumor_Sample_Barcode %in% old_samples)

# compile maf
young_maf <- read.maf(maf = young_maf, clinicalData = young_clinical)
old_maf <- read.maf(maf = old_maf, clinicalData = old_clinical)

# compare mafs --> forest plot
compare_maf <- mafCompare(m1 = young_maf, m2 = old_maf, m1Name = "Young", m2Name = "Old", minMut = 5)

# co oncoplot
gns <- c("PTEN", "PIK3CA", "TP53", "CTCF", "KRAS")
tiff("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Figures/updated_mutation/ucec_co_oncoplot.tiff", res = 320,
     width = 6, height = 3, units = "in")
coOncoplot(m1 = young_maf, m2 = old_maf, genes = gns, m1Name = "Young", m2Name = "Old",sepwd_genes2 = 0,
           sepwd_genes1 = 0)
dev.off()

# co lollipop
tiff("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Figures/updated_mutation/ucec_co_lollipop.tiff", res = 320,
     width = 10, height = 4, units = "in")
lollipopPlot2(m1 = young_maf, m2 = old_maf, gene = "PTEN", m1_name = "Young", m2_name = "Old",
              AACol1 = "HGVSp", AACol2 = "HGVSp", legendTxtSize = 2, pointSize = 2, domainLabelSize = 1.5)
dev.off()


AACol2 = # some pathways
pdf("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Figures/updated_mutation/ucec_mut_pathways.pdf")
maftools::PlotOncogenicPathways(young_maf, "RTK-RAS")
legend("topright", "Young", bty="n", cex = 2)
maftools::PlotOncogenicPathways(old_maf, "RTK-RAS")
legend("topright", "Old", bty="n", cex = 2)
maftools::PlotOncogenicPathways(young_maf, "PI3K")
legend("topright", "Young", bty="n", cex = 2)
maftools::PlotOncogenicPathways(old_maf, "PI3K")
legend("topright", "Old", bty="n", cex = 2)
dev.off()


# compare boxplots
young <- data.frame(young_maf@variant.classification.summary)
young$Age <- "Young"
old <- data.frame(old_maf@variant.classification.summary)
old$Age <- "Old"
combined <- rbind(young, old)
colnames(combined) <- gsub('\\_',' ',colnames(combined))
combined <- combined[, -c(1,11)]
combined <- melt(combined)
combined$Age <- factor(combined$Age, levels = c("Young", "Old"))
compareTest <- combined %>%
  dplyr::group_by(variable) %>%
  rstatix::wilcox_test(data = ., formula = value ~ Age)
compareTest$p.adj <- p.adjust(compareTest$p, method = "fdr")
compareTest$fdrstar <- stars.pval(compareTest$p.adj)

ggboxplot(combined, x = 'Age', y = 'value', fill = 'Age', palette = 'jco', ylab = "Frequency",
          title = "UCEC SNPs") + 
  facet_wrap(~variable, labeller = label_wrap_gen(width=15), scales = "free_y") + scale_y_log10() + 
  stat_pvalue_manual(data = compareTest, y.position = c(2.5,3.2,1.5,3,1.5), 
                                       label = "fdrstar", hide.ns = T, size = 15, vjust = 1, bracket.size = 0) +
  theme_pubr(base_size = 35) + theme(axis.title.x = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     strip.text = element_text(size = 25))
ggsave("/Users/Yajas/Documents/Elemento/AgeTCGA-master/Final/Figures/updated_mutation/ucec_variants.tiff",
       width = 12, height = 12)




## modify maftools scripts
mod.coOncoplot = function(m1, m2, genes = NULL, m1Name = NULL, m2Name = NULL,
                      clinicalFeatures1 = NULL, clinicalFeatures2 = NULL,
                      annotationColor1 = NULL, annotationColor2 = NULL, annotationFontSize = 1.2,
                      sortByAnnotation1 = FALSE, sortByAnnotation2 = FALSE,
                      additionalFeature1 = NULL, additionalFeaturePch1 = 20, additionalFeatureCol1 = "white", additionalFeatureCex1 = 0.9,
                      additionalFeature2 = NULL, additionalFeaturePch2 = 20, additionalFeatureCol2 = "white", additionalFeatureCex2 = 0.9,
                      sepwd_genes1 = 0.5, sepwd_samples1 = 0.5, sepwd_genes2 = 0.5, sepwd_samples2 = 0.5,
                      colors = NULL, removeNonMutated = TRUE,
                      geneNamefont = 0.8, showSampleNames = FALSE, SampleNamefont = 1,
                      legendFontSize = 1.2, titleFontSize = 1.5, keepGeneOrder=FALSE,
                      bgCol = "#CCCCCC", borderCol = "white"){
  
  if(is.null(genes)){
    m1.genes = getGeneSummary(m1)[1:5]
    m2.genes = getGeneSummary(m2)[1:5]
    mdt = merge(m1.genes[,.(Hugo_Symbol, MutatedSamples)], m2.genes[,.(Hugo_Symbol, MutatedSamples)], by = 'Hugo_Symbol', all = TRUE)
    mdt$MutatedSamples.x[is.na(mdt$MutatedSamples.x)] = 0
    mdt$MutatedSamples.y[is.na(mdt$MutatedSamples.y)] = 0
    mdt$max = apply(mdt[,.(MutatedSamples.x, MutatedSamples.y)], 1, max)
    mdt = mdt[order(max, decreasing = TRUE)]
    
    genes = mdt[,Hugo_Symbol]
  }
  
  m1.sampleSize = m1@summary[3, summary]
  m2.sampleSize = m2@summary[3, summary]
  
  
  if(is.null(m1Name)){
    m1Name = 'M1'
  }
  
  m1Name = paste(m1Name, ' (N = ' , m1.sampleSize, ')',sep = '')
  
  if(is.null(m2Name)){
    m2Name = 'M2'
  }
  
  m2Name = paste(m2Name, ' (N = ' , m2.sampleSize, ')',sep = '')
  
  if(is.null(colors)){
    vc_col = get_vcColors()
  }else{
    vc_col = colors
  }
  
  
  m12_annotation_colors = NULL
  if(!is.null(clinicalFeatures1) & !is.null(clinicalFeatures2)){
    if(is.null(annotationColor1) & is.null(annotationColor2)){
      m12_annotation_colors = get_m12_annotation_colors(a1 = m1, a1_cf = clinicalFeatures1,
                                                        a2 = m2, a2_cf = clinicalFeatures2)
      annotationColor1 = m12_annotation_colors
      annotationColor2 = m12_annotation_colors
    }
  }
  
  #Get matrix dimensions and legends to adjust plot graphics::layout
  nm1 = print_mat(maf = m1, genes = genes, removeNonMutated = removeNonMutated,
                  test = TRUE, colors = colors)
  nm1_ncol = ncol(nm1[[1]])
  nm1_vc_cols = nm1[[2]]
  
  nm2 = print_mat(maf = m2, genes = genes, removeNonMutated = removeNonMutated,
                  test = TRUE, colors = colors)
  nm2_ncol = ncol(nm2[[1]])
  nm2_vc_cols = nm2[[2]]
  
  if(!is.null(clinicalFeatures1) || !is.null(clinicalFeatures2)){
    mat_lo = mat_lo = matrix(data = c(1,3,5,2,4,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
    mat_lo = graphics::layout(mat = mat_lo,
                              widths = c(6 * (nm1_ncol/nm2_ncol), 1.5, 6), heights = c(12, 3, 4))
  }else{
    mat_lo = matrix(data = c(1,2,3,4,4,4), nrow = 2, ncol = 3, byrow = TRUE)
    mat_lo = graphics::layout(mat = mat_lo,
                              widths = c(6 * (nm1_ncol/nm2_ncol), 1, 6), heights = c(12, 4))
  }
  
  
  m1_legend = print_mat(maf = m1, genes = genes, removeNonMutated = removeNonMutated,
                        clinicalFeatures = clinicalFeatures1, colors = colors,
                        annotationColor = annotationColor1, barcode_size = SampleNamefont,
                        sortByAnnotation = sortByAnnotation1, fontSize = geneNamefont,
                        title = m1Name, title_size = titleFontSize,
                        showBarcodes = showSampleNames, bgCol = bgCol, borderCol = borderCol,
                        additionalFeature = additionalFeature1, additionalFeaturePch = additionalFeaturePch1,
                        additionalFeatureCex = additionalFeatureCex1, additionalFeatureCol = additionalFeatureCol1, sepwd_genes = sepwd_genes1, sepwd_samples = sepwd_samples1)
  
  
  if(is.null(clinicalFeatures1) & !is.null(clinicalFeatures2)){
    plot.new()
  }
  
  if(showSampleNames){
    par(mar = c(5, 0, 3, 0))
  }else{
    par(mar = c(1, 0, 3, 0))
  }
  
  nm = matrix(data = 1, nrow = 1, ncol = length(genes))
  image(x = 1:nrow(nm), y = 1:ncol(nm), z = nm, axes = FALSE, xaxt="n", yaxt="n",
        xlab="", ylab="", col = "white")
  mtext(text = rev(genes), side = 2, adj = 0.5, at = 1:ncol(nm),
        font = 3, line = -2, cex = geneNamefont, las = 2)
  
  if(!is.null(clinicalFeatures1) || !is.null(clinicalFeatures2)){
    plot.new()
  }
  
  m2_legend = print_mat(maf = m2, genes = genes, removeNonMutated = removeNonMutated,
                        clinicalFeatures = clinicalFeatures2, colors = colors,
                        annotationColor = annotationColor2, barcode_size = SampleNamefont,
                        sortByAnnotation = sortByAnnotation2, fontSize = geneNamefont,
                        title = m2Name, title_size = titleFontSize, plot2 = TRUE,
                        showBarcodes = showSampleNames, bgCol = bgCol, borderCol = borderCol,
                        additionalFeature = additionalFeature2, additionalFeaturePch = additionalFeaturePch2,
                        additionalFeatureCex = additionalFeatureCex2, additionalFeatureCol = additionalFeatureCol2,
                        sepwd_genes = sepwd_genes2, sepwd_samples = sepwd_samples2)
  
  if(!is.null(clinicalFeatures1) & is.null(clinicalFeatures2)){
    plot.new()
  }
  
  vc_legend = unique(c(names(nm1_vc_cols), names(nm2_vc_cols)))
  vc_legend = vc_col[vc_legend]
  vc_legend = vc_legend[!is.na(vc_legend)]
  
  
  if(is.null(m12_annotation_colors)){
    anno_legend = c(m1_legend, m2_legend)
  }else{
    anno_legend = m12_annotation_colors
  }
  
  par(mar = c(2, 2, 0, 0), xpd = TRUE)
  
  vc_pch = rep(15, length(vc_legend))
  if(!is.null(additionalFeature1)){
    vc_legend = c(vc_legend, "gray70")
    names(vc_legend)[length(vc_legend)] = paste(additionalFeature1, collapse = ":")
    vc_pch = c(vc_pch, additionalFeaturePch1)
  }
  
  if(!is.null(additionalFeature2)){
    vc_legend = c(vc_legend, "gray70")
    names(vc_legend)[length(vc_legend)] = paste(additionalFeature2, collapse = ":")
    vc_pch = c(vc_pch, additionalFeaturePch2)
  }
  
  plot(NULL,ylab='',xlab='', xlim=0:1, ylim=0:1, axes = FALSE)
  lep = legend("topleft", legend = names(vc_legend),
               col = vc_legend, border = NA, bty = "n",
               ncol= 2, pch = vc_pch, xpd = TRUE, xjust = 0, yjust = 0, cex = legendFontSize)
  
  x_axp = 0+lep$rect$w
  
  if(!is.null(anno_legend)){
    
    for(i in 1:length(anno_legend)){
      #x = unique(annotation[,i])
      x = anno_legend[[i]]
      
      if(length(x) <= 4){
        n_col = 1
      }else{
        n_col = (length(x) %/% 4)+1
      }
      
      lep = legend(x = x_axp, y = 1, legend = names(x),
                   col = x, border = NA,
                   ncol= n_col, pch = 15, xpd = TRUE, xjust = 0, bty = "n",
                   cex = annotationFontSize, title = names(anno_legend)[i], title.adj = 0)
      x_axp = x_axp + lep$rect$w
      
    }
  }
  #title()
}
get_vcColors = function(alpha = 1){
  col = c(RColorBrewer::brewer.pal(11, name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue', '#7b7060')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                         'RNA','Splice_Site','Intron','Frame_Shift_Ins','In_Frame_Del','ITD','In_Frame_Ins',
                         'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del', 'Complex_Event')
  col
}
