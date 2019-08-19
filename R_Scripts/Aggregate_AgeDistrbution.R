##########################
#Task: Go through all the clinical files and estimate the age distribution
# Notes/Updates
##########################


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
    Medi_Y = nrow(clin_data_df[which(clin_data_df$AGE < Median_Age),])
    Medi_O = nrow(clin_data_df[which(clin_data_df$AGE > Median_Age),])
    Mean_Y = nrow(clin_data_df[which(clin_data_df$AGE < Mean_Age),])
    Mean_O = nrow(clin_data_df[which(clin_data_df$AGE > Mean_Age),])
    TotalNumberCases =
    
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



# Append to dataframe above the median_younger/median_older columns 
#  this would be samples < median age , sample > median age
# same this for the mean_younger, mean_older


#you should end up with 2 object here a list object for the clinical table and dataframe summary to determine the group distributions
