
########################################################################
####          Subset common samples between both datasets           ####
########################################################################

library(TCGAutils)
filtered_counts <- RNAseq_counts
filtered_counts <- filtered_counts[ ,(TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Solid Tumor")]    ## only look at primary solid tumors
which(duplicated(TCGAbiospec(colnames(filtered_counts))$submitter_id))    ## there are duplicate case IDs. remove all of them
filtered_counts <- filtered_counts[,-c(which(duplicated(TCGAbiospec(colnames(filtered_counts))$submitter_id)))]    ## removing duplicates
colnames(filtered_counts) <- substr(colnames(filtered_counts), start = 1, stop = 15)    ## making sure both datasets have the same case id format
length(which(colnames(filtered_counts) %in% Age_OnyCT.df.Y$CASE_ID))    ## 8878 subjects for analysis
filtered_age <- Age_OnyCT.df.Y[which(Age_OnyCT.df.Y$CASE_ID %in% colnames(filtered_counts)),]




########################################################################


########################################################################
####                Subtype Specific Age Based (Q1/Q3)              ####
####                    Differential Expression                     ####
########################################################################

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(biomaRt)

cancer_types <- levels(factor(Age_OnyCT.df$TCGA_CancerCode))    ## unique 3 letter codes


#### loop through all cancer types to generate DE + volcano plots
for (i in length(cancer_types)) {
  
  #### create a subset of the age and rna seq datasets
  #### select the current cancer type subjects which have which are either <Q1 or>Q3 wrt age
  age_subset <- filtered_age[which(filtered_age$TCGA_CancerCode == cancer_types[i] & (filtered_age$Subtype_First_Quartile == "Young_Q1" | filtered_age$Subtype_Last_Quartile == "Old_Q3")),]
  rna_subset <- filtered_counts[,which(colnames(filtered_counts) %in% age_subset$CASE_ID)]
  
  #### generate a data frame with case ID and designating subsetted samples to either "Young_Q1" or "Old_Q3"
  input <- as.data.frame(cbind(age_subset$CASE_ID, age_subset$Subtype_First_Quartile, age_subset$Subtype_Last_Quartile))
  
  #### loop through the data frame so as to combine Q1- and Q3+ into a single column (V4)
  for (i in 1:nrow(input)) {
    input[i,4] <- ifelse(input$V2[i] == "Young_Q1", "Young_Q1", (ifelse(input$V3[i] == "Old_Q3", "Old_Q3", NA)))
  }
  input <- na.omit(input)
  input <- input[,-c(2,3)]
  
  #### select the subjects present in the dataframe created above
  rna_subset <- filtered_counts[,which(colnames(filtered_counts) %in% input$V1)]
  
  
  
  #### Begin DE using edgeR
  dgList <- DGEList(counts=rna_subset, genes=rownames(rna_subset))  ## create a DGEList object
  keep <- which(rowSums(rna_subset) >= 2)                           ## identify rows which have >2 counts/row
  dgList <- dgList[keep,]                                           ## keep only the rows which have >2 counts/row
  designMat <- model.matrix(~0 + factor(input$V4))                  ## remember that input$V4 is a vector with Q1-/Q3+ information
  colnames(designMat) <- c("Old_Q3", "Young_Q1")                    ## change design matrix column names
  dgList <- estimateDisp(dgList, designMat, robust=TRUE)            ## estimate dispersions. can be viewed using plotBCV()
  fit <- glmQLFit(dgList, designMat, robust=TRUE)                   ## model fitting 1
  contr.mat <- makeContrasts(Old_Q3-Young_Q1, levels=designMat)     ## create a contrast matrix
  res <- glmQLFTest(fit, contrast=contr.mat)                        ## model fitting 2
  tt <- data.frame(topTags(res,n = nrow(rna_suset)))                ## toptable with DE data
  
  
  ## create a new dataframe with significance annotations. used in the volcano plot
  mutateddf1 = mutate(tt, Significance=ifelse(tt$FDR<0.05, "FDR<0.05", "Not Sig"))
  
  
  
  #### biomaRt. it would save a lot of time if this was done outside the loop
  mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")  ## convert ensembl to gene name
  #### dat is a 2 column data frame with ensembl and gene name columns. Use this to match the ensembl DE results
  dat = getBM(
    values = mutateddf1$genes,
    filters = c("ensembl_gene_id"),
    attributes = c("ensembl_gene_id", "external_gene_name"),
    mart = mart
  )
  colnames(dat)[1] <- "genes"
  
  
  ## going back to the new data frame with significance annotations
  mutateddf1 <- merge(x = mutateddf1, y = dat, by = 'genes')        ## add the gene id column
  mutateddf1 <- mutateddf1[order(mutateddf1$FDR),]                  ## reorder by FDR
  
  
  ## Volcano plot
  volc = ggplot(mutateddf1, aes(logFC, -log10(FDR))) +
    geom_point(aes(col=Significance)) +
    scale_color_manual(values=c("red", "black")) + 
    theme(legend.title = element_blank())+
    xlab(bquote(''~log[2]~'Fold Change'))+
    ylab(bquote(''~-log[10]~~italic(FDR)~ ''))+
    theme(legend.position = "none")+
    ggtitle(paste(cancer_types[i],"    Old_Q3-Young_Q1"))
  volc=volc+geom_text_repel(data=head(mutateddf1, 25), aes(label=external_gene_name))
  
  ggsave(filename =  paste0(cancer_types[i]," Quartile DE.jpg"),plot = volc,dpi = 320)
}


## talk about getting the loop to work
