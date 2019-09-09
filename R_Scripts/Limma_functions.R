Create_Design<-function(input,condition,comparison){
  library("limma")
 # print(condition)
  input<- as.data.frame(input)
  colnames(input)<-c("Samples","target")
  print(head(input))
  input$target <- as.matrix(input$target)
    design <<- model.matrix(~0+ unlist(input$target))
  colnames(design) <<- condition
  print(head(design))
  print(as.character(comparison))
    contmatrix <<- makeContrasts(as.character(comparison),levels=design)
}
CallLimma<-function(Data,design,contrastMatrix){ 
    suppressMessages(library(limma))
  suppressMessages(library(edgeR))
  print(dim(Data))
  Data<-Data[which(rowSums(cpm(Data)> 1) >= 2) ,]
  print(dim(Data))
    Data.voom<- voom(Data,plot=FALSE)
    Voom_output <<- Data.voom
    fit<- lmFit(Data.voom,design)
    fit2<<- contrasts.fit(fit,contrastMatrix)
    fit2<<-eBayes(fit2)
    return(topTable(fit2,n=Inf))
    
}


### To run DE Code ###
input <-Urine_all_batches_Kept[which(Urine_all_batches_Kept$Diagnosis %in% c("NoRejection","AMR")),]  # Annotation
input <- as.data.frame(input[,c(1,6)])
condition <- c("AMR","NoRejection")
comparison <- paste(condition[1],condition[2],sep="-")
Create_Design(input,condition,comparison )
Count_input <- Combined_Batches_UrineRNASeq[,as.matrix(input$Sample_Name)]
DE_fil_AMR_vs_Normal <- CallLimma(Count_input,design,contmatrix)
DE_fil_AMR_vs_Normal <-merge(All_Urine_Samples_RNASeq_done[,1,drop=FALSE],DE_fil_AMR_vs_Normal,by='row.names')
write.xlsx2(DE_fil_AMR_vs_Normal[which(DE_fil_AMR_vs_Normal$adj.P.Val < 0.1 & DE_fil_AMR_vs_Normal$logFC > 0),],
            'Diff_Exp_AMR_vs_Normal_June_2018.xlsx',sheetName = 'Up AMR vs Normal',row.names = FALSE)
write.xlsx2(DE_fil_AMR_vs_Normal[which(DE_fil_AMR_vs_Normal$adj.P.Val < 0.1  & DE_fil_AMR_vs_Normal$logFC < 0),],
            'Diff_Exp_AMR_vs_Normal_June_2018.xlsx',sheetName = 'Down AMR vs Normal',append = 1,,row.names = FALSE)