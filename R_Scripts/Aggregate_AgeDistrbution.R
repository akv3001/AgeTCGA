##########################
#Task: Go through all the clinical files and estimate the age distribution
# Notes/Updates
##########################

library(data.table)
library(ggpubr)
# Write a loop to iterate through the clinical files 
# Create a list object with each cancer type as adataframe object

path = "DATA/cBioportal_Survival_Data"  ## set file path 

file.names <- list.files(path, pattern =".txt",full.names = TRUE) ##assign name to path (looking for txt files withi directory)

All_CancerDataClinicalData <- list()    ##open list
Summarized_Ages.AllCancers =  data.frame(CancerType=character(0),Median=numeric(0),Mean=numeric(0)
                                         ,Max=numeric(0),
                                 Min=numeric(0),Median_Younger=numeric(0),Median_Old=numeric(0),
                                 Mean_Younger=numeric(0), Mean_Older=numeric(0))


for (i in 1:length(file.names)) ## for object(i) in 1 through length of file.names. 
{  
    print(i)
    clin_data_df <- read.table(file.names[[i]],sep="\t",stringsAsFactors = FALSE,strip.white = TRUE,fill = TRUE,header = TRUE,blank.lines.skip = TRUE)  ##create dataframe of cBio files for all cancer files
    CancerType = unlist(strsplit(basename(file.names[[i]]),split = "_"))[[3]]
    print(CancerType)
  
    ##save cancer typeto var cancerType 
    clin_data_df$AGE = as.numeric(clin_data_df$AGE)
    Median_Age=median(clin_data_df$AGE,na.rm = TRUE)
    Mean_Age = mean(clin_data_df$AGE,na.rm = TRUE)
    Min_Age = min(clin_data_df$AGE,na.rm = TRUE)
    Max_Age = max(clin_data_df$AGE,na.rm = TRUE)
    
    # Append to dataframe above the median_younger/median_older columns 
    #  this would be samples < median age , sample > median age
    # same this for the mean_younger, mean_older
    Medi_Y = nrow(clin_data_df[which(clin_data_df$AGE < Median_Age),])
    Medi_O = nrow(clin_data_df[which(clin_data_df$AGE > Median_Age),])
    Mean_Y = nrow(clin_data_df[which(clin_data_df$AGE < Mean_Age),])
    Mean_O = nrow(clin_data_df[which(clin_data_df$AGE > Mean_Age),])
    TotalNumberCases =nrow(clin_data_df)
    
    clin_data_df$Median_Classification = ifelse(clin_data_df$AGE > Median_Age ,"Old_Median","Young_Median")
    clin_data_df$Mean_Classification = ifelse(clin_data_df$AGE > Mean_Age , "Old_Mean","Young_Mean")
    
    
    if(Median_Age > 60){
      clin_data_df$CancerAgeMedianCategory = "Later-Life"
    }
    else {
      clin_data_df$CancerAgeMedianCategory = "Younger-Life"
      
    }
    # Compute the median and mean for each cancer type , the mix , max 
    # also use a loop to do this by going through the list structure
    # this should create a new dataframe summarized ( you did this manually in the excel originally but please execute same in code below)
    
    Summarized_Ages.AllCancers = rbind(Summarized_Ages.AllCancers,
                    data.frame(CancerType=CancerType,Median=Median_Age,Mean = Mean_Age,
                               Max=Max_Age,Min=Min_Age,
                               Median_Younger=Medi_Y,Median_Old=Medi_O,
               Mean_Younger=Mean_Y, Mean_Older=Mean_O)
    )

    All_CancerDataClinicalData[[CancerType]] <- clin_data_df ## save to new list 
}

##print list
print(All_CancerData) 
print(AgeDistribution_Agegroup) 

### Merge Ages across all cancers ##
## Mostly needed for visualization ##
Age_OnyCT= list()
for ( EachCT in 1:length(All_CancerDataClinicalData)){
  name_CT = names(All_CancerDataClinicalData)[[EachCT]]
  print(name_CT)
  selCT = All_CancerDataClinicalData[[name_CT]]
  selCT = selCT[,c("CASE_ID","AGE","CANCER_TYPE","Median_Classification","Mean_Classification","CancerAgeMedianCategory")]
  selCT$TCGA_CancerCode = name_CT
  Age_OnyCT[[name_CT]] = selCT
}

Age_OnyCT.df= rbindlist(Age_OnyCT)
Age_OnyCT.df=na.omit(Age_OnyCT.df)

ggviolin(Age_OnyCT.df,x = 'TCGA_CancerCode',y= 'AGE' ,fill = 'TCGA_CancerCode',
          order = Summarized_Ages.AllCancers$CancerType,add = 'jitter')+ geom_hline(yintercept = 60,color="red")


Summarized_Ages.AllCancers = Summarized_Ages.AllCancers[order(Summarized_Ages.AllCancers$Median,decreasing = TRUE),]
ggboxplot(Age_OnyCT.df,x = 'TCGA_CancerCode',y= 'AGE' ,
          fill = 'CancerAgeMedianCategory',palette = 'jco',
         order = Summarized_Ages.AllCancers$CancerType)+
         #add = 'jitter' 
  geom_hline(yintercept = 60,color="red")


# Comparing Mean Age across all cancers
pdf('Comparison_MeanAge_AcrossAllCancersBoxplot.pdf')
compare_stat = list(c("Young_Mean","Old_Mean"))
print(ggboxplot(Age_OnyCT.df,x = 'Mean_Classification',y= 'AGE' ,
          fill = 'Mean_Classification',palette = 'jco')+
  stat_compare_means(comparisons = compare_stat,method = "wilcox"))
dev.off()

pdf('Comparison_MedianAge_AcrossAllCancersBoxplot.pdf')
compare_stat = list(c("Young_Median","Old_Median"))
print(ggboxplot(Age_OnyCT.df,x = 'Median_Classification',y= 'AGE' ,
                fill = 'Median_Classification',palette = 'jco')+
        stat_compare_means(comparisons = compare_stat,method = "wilcox"))
dev.off()

## create new excel file from df Summarized_Ages.AllCancers 
library("openxlsx")
write.xlsx(Summarized_Ages.AllCancers, file = "Age_CancerType_Summary.xlsx", colNames = TRUE, borders = "columns")

## Choose Top 5 younger age Cancer Types by (median age) create new DF ##create new pdf & excel file with Top 5 young age CT
young_CT <- head(Summarized_Ages.AllCancers[order(Summarized_Ages.AllCancers$Median, decreasing= F),], n = 5)
pdf("younger5_CT.pdf")
write.xlsx(young_CT, file = "younger5_CT.xlsx", colNames = TRUE, borders = "columns")
dev.off()

## Choose Top 5 older age Cancer Types by (median age) create new DF ##create new pdf & excel file with Top 5 older age CT
older_CT <- head(Summarized_Ages.AllCancers[order(Summarized_Ages.AllCancers$Median, decreasing= T),], n = 5)
pdf("younger5_CT.pdf")
write.xlsx(older_CT, file = "older5_CT.xlsx", colNames = TRUE, borders = "columns")
dev.off()

############################Test plot################ 
mdold <- Summarized_Ages.AllCancers$Median_Old 
mdyoung <- Summarized_Ages.AllCancers$Median_Younger
qplot(mdold, mdyoung, data=Summarized_Ages.AllCancers, color= CancerType) 


#you should end up with 2 object here a list object for the
#clinical table and dataframe summary to determine the group distributions

