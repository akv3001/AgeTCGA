
## classify subjects based on their age and cancer subtypes
## Classifiers are either > or < these: mean, median,Q1,Q3
## classification occurs within each cancer subtype

########################################################################
####                      Subtype Classification                    ####
########################################################################


cancer_types <- levels(factor(Age_OnyCT.df$TCGA_CancerCode))    ## unique 3 letter codes

#### setup a data frame of summary stats for all cancers
descriptive_stats <- matrix(nrow = length(cancer_types), ncol = 6)
rownames(descriptive_stats) <- cancer_types
colnames(descriptive_stats) <- names(summary(Age_OnyCT.df$AGE))
descriptive_stats <- as.data.frame(descriptive_stats)

#### loop through cancer types to populate the summary stats data frame
for (i in 1:length(cancer_types)) {
  row_index <- which(Age_OnyCT.df$TCGA_CancerCode == cancer_types[i])
  descriptive_stats[i,] <- summary(Age_OnyCT.df[row_index,]$AGE)
}
descriptive_stats    ## inspect

#### setup a dataframe for classifying subjects within cancer type
#### cbind this to Age_OnyCT.df once you fill it up
new_classification <- matrix(ncol = 4)
colnames(new_classification) <- c("Subtype_Mean_Classification", "Subtype_Median_Classification", "Subtype_First_Quartile", "Subtype_Last_Quartile")
new_classification <- as.data.frame(new_classification)

#### loop to populate the classification table
#### quartile limits are included
#### > mean/median
for (i in 1:nrow(Age_OnyCT.df)) {
  for (j in 1:nrow(descriptive_stats)) {
    if (Age_OnyCT.df[i,]$TCGA_CancerCode == rownames(descriptive_stats)[j]) {
      new_classification[i,]$Subtype_Median_Classification = ifelse(Age_OnyCT.df[i,]$AGE > descriptive_stats[j,]$Median,
                                                                    "Old_Median",
                                                                    "Young_Median")
      new_classification[i,]$Subtype_Mean_Classification = ifelse(Age_OnyCT.df[i,]$AGE > descriptive_stats[j,]$Mean,
                                                                    "Old_Mean",
                                                                    "Young_Mean")
      new_classification[i,]$Subtype_First_Quartile = ifelse(Age_OnyCT.df[i,]$AGE > descriptive_stats[j,]$`1st Qu.`,
                                                                    "Old_Q1",
                                                                    "Young_Q1")
      new_classification[i,]$Subtype_Last_Quartile = ifelse(Age_OnyCT.df[i,]$AGE >= descriptive_stats[j,]$`3rd Qu.`,
                                                                    "Old_Q3",
                                                                    "Young_Q3")
    }
  }
}
head(new_classification)    ## inspect

Age_OnyCT.df.Y <- cbind(Age_OnyCT.df,new_classification)    ## create a single DF with all variables
head(Age_OnyCT.df.Y)    ## inspect





## classify subjects based on their age
## Classifiers are either > or < these: Q1,Q3
## classification occurs globally
########################################################################
####                          PanCan Q1/Q3                          ####
########################################################################

#### Similarly, add Q1 and Q3 for the overall classification (PanCan)
#### setup a dataframe for classifying subjects within cancer type
#### cbind this to Age_OnyCT.df once you fill it up
new_classification <- matrix(ncol = 2)
colnames(new_classification) <- c("PanCan_First_Quartile", "PanCan_Last_Quartile")
new_classification <- as.data.frame(new_classification)

#### loop to populate the classification table
#### quartile limits are included
for (i in 1:nrow(Age_OnyCT.df)) {
  new_classification[i,]$PanCan_First_Quartile = ifelse(Age_OnyCT.df[i,]$AGE > summary(Age_OnyCT.df$AGE)[2],
                                                         "Old_Q1",
                                                         "Young_Q1")
  new_classification[i,]$PanCan_Last_Quartile = ifelse(Age_OnyCT.df[i,]$AGE >= summary(Age_OnyCT.df$AGE)[5],
                                                        "Old_Q3",
                                                        "Young_Q3")
}
head(new_classification)  ## inspect

Age_OnyCT.df.Y <- cbind(Age_OnyCT.df.Y, new_classification)    ## create a single DF with all variables
head(Age_OnyCT.df.Y)    ## inspect
colnames(Age_OnyCT.df.Y)
##
library(ggpubr)
ggviolin(data = Age_OnyCT.df.Y, x = 'TCGA_CancerCode', y = 'AGE', add = 'jitter', 
         color = c('PanCan_Last_Quartile'))



## lines 53 and 90 show NA for the second rowname - why?



