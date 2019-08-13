##########################
#Task: Go through all the clinical files and estimate the age distribution
# Notes/Updates


##########################


# Write a loop to iterate through the clinical files 
# Create a list object with each cancer type as adataframe object

path = "DATA/cBioportal_test"  ## set file path 

file.names <- list.files(path, pattern =".txt",full.names = TRUE) ##assign name to path (looking for txt files withi directory)


## for object(i) in 1 through length of file.names

for (i in 1:length(file.names)) ## for object(i) in 1 through length of file.names
  
{
  
  print(file.names[[i]])
  
  
}


# Compute the median and mean for each cancer type , the mix , max 
# also use a loop to do this by going through the list structure
# this should create a new dataframe summarized ( you did this manually in the excel originally but please execute same in code below)



# Append to dataframe above the median_younger/median_older columns 
#  this would be samples < median age , sample > median age
# same this for the mean_younger, mean_older



#you should end up with 2 object here a list object for the clinical table and dataframe summary to determine the group distributions
