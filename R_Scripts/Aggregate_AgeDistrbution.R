##########################
#Task: Go through all the clinical files and estimate the age distribution
# Notes/Updates
##########################


# Write a loop to iterate through the clinical files 
# Create a list object with each cancer type as adataframe object

path = "DATA/cBioportal_test"  ## set file path 

file.names <- list.files(path, pattern =".txt",full.names = TRUE) ##assign name to path (looking for txt files withi directory)

All_CancerData <- list()    ##open list

for (i in 1:length(file.names)) ## for object(i) in 1 through length of file.names. 
{  
  
  clin_data_df <- read.delim(file.names[[i]], stringsAsFactors = F)  ##create dataframe of cBio files for all cancer files
  
  
  for ( type in 1:length(unique(clin_data_df$CANCER_TYPE))){  ##loop through the number cancer type
    
    CancerType <- (clin_data_df$CANCER_TYPE)[type]   ##save cancer typeto var cancerType 
    SelectCases <- na.omit(clin_data_df$AGE[which(clin_data_df$CANCER_TYPE == CancerType)]) ## extracting age for the cancerType
    
    
    # Compute the median and mean for each cancer type , the mix , max 
    # also use a loop to do this by going through the list structure
    # this should create a new dataframe summarized ( you did this manually in the excel originally but please execute same in code below)
    
    
    ##calculating range/quantile/median/mean/max/min/greater > median/ less than <median 
    CancerData <- c(range=range(SelectCases), 
                    median=median(SelectCases),
                    mean=mean(SelectCases), min = SelectCases[sapply(which.min(SelectCases),min)],
                    max=SelectCases[sapply(which.max(SelectCases),max)],
                    OlderMedian=sum(SelectCases>median(SelectCases)),
                    YoungerMedian=sum(SelectCases<median(SelectCases)),
                    quantile=quantile(SelectCases))
    
    All_CancerData[[CancerType]] <- CancerData ## save to new list 
    
  }
}

print(All_CancerData) ##print list 

# Append to dataframe above the median_younger/median_older columns 
#  this would be samples < median age , sample > median age
# same this for the mean_younger, mean_older


#you should end up with 2 object here a list object for the clinical table and dataframe summary to determine the group distributions
