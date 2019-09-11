####################
# Task :Use TCGA biolinks to download htseq counts data

#####################
library(TCGAbiolinks)


All_queries.list <- list()

for (i in 1:length(file.names)) ## for object(i) in 1 through length of file.names. 
{  
  print(i)
  CancerType = unlist(strsplit(basename(file.names[[i]]),split = "_"))[[3]]
  
  print(CancerType)
  
  CancerType=paste('TCGA-',toupper(CancerType),sep="")
  query <- GDCquery(project=CancerType, 
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type   = "HTSeq - Counts",
                    experimental.strategy = "RNA-Seq",)
  GDCdownload(query)
  
  All_queries.list[[CancerType]] <- query
  
}
print('Done all')

save(All_queries.list,'All_queries_htseqCounts_list.rds')
####### Download Matrix ###
TCGA_hg38_Cancers_HTSeqCount <- list()
for (CancerType in names(All_queries.list)) ## for object(i) in 1 through length of file.names. 
{  
  print(CancerType)

  TCGA_hg38_Cancers_HTSeqCount[[CancerType]] <-  GDCprepare(All_queries.list[[CancerType]],summarizedExperiment = F)
 # colnames(TCGA_hg38_Cancers_HTSeqCount[[CancerType]]) <- gsub("normalized_count_","",colnames(TCGA_all_Cancers_hg19_Normalized_Expression[[CancerType]]))
}

saveRDS(TCGA_hg38_Cancers_HTSeqCount,'TCGA_hg38_Cancers_HTSeqCount.rds')




