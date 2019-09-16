
## classify subjects based on their age and cancer subtypes
## Classifiers are either > or < these: mean, median,Q1,Q3
## classification occurs within each cancer subtype

########################################################################
####                      Subtype Classification                    ####
########################################################################

cancer_types <- levels(factor(Age_OnyCT.df$TCGA_CancerCode))    ## unique 3 letter codes

#### setup a data frame of summary stats for all cancers
descriptive_stats <- matrix(nrow = length(cancer_types), ncol = 8)
rownames(descriptive_stats) <- cancer_types
colnames(descriptive_stats) <- c(names(summary(Age_OnyCT.df$AGE)), "T1", "T3")
descriptive_stats <- as.data.frame(descriptive_stats)

#### loop through cancer types to populate the summary stats data frame
for (i in 1:length(cancer_types)) {
  row_index <- which(Age_OnyCT.df$TCGA_CancerCode == cancer_types[i])
  descriptive_stats[i,1:6] <- summary(Age_OnyCT.df[row_index,]$AGE)
  descriptive_stats[i, 7] <- quantile(Age_OnyCT.df[row_index,]$AGE, prob = 0.33)
  descriptive_stats[i, 8] <- quantile(Age_OnyCT.df[row_index,]$AGE, prob = 0.66)
}
descriptive_stats    ## inspect

#### setup a dataframe for classifying subjects within cancer type
#### cbind this to Age_OnyCT.df once you fill it up
new_classification <- matrix(ncol = 6)
colnames(new_classification) <- c("Subtype_Mean_Classification", "Subtype_Median_Classification", "Subtype_First_Quartile", "Subtype_Last_Quartile", "Subtype_First_Tertile", "Subtype_Last_Tertile")
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
      new_classification[i,]$Subtype_First_Tertile = ifelse(Age_OnyCT.df[i,]$AGE <= descriptive_stats[j,]$T1,
                                                            "Young_T1",
                                                            "Old_T1")
      new_classification[i,]$Subtype_Last_Tertile = ifelse(Age_OnyCT.df[i,]$AGE >= descriptive_stats[j,]$T3,
                                                            "Old_T3",
                                                            "Young_T3")
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
new_classification <- matrix(ncol = 4)
colnames(new_classification) <- c("PanCan_First_Quartile", "PanCan_Last_Quartile", "PanCan_First_Tertile", "PanCan_Last_Tertile")
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
  new_classification[i,]$PanCan_First_Tertile = ifelse(Age_OnyCT.df[i,]$AGE <= quantile(Age_OnyCT.df$AGE, prob = 0.33),
                                                         "Young_T1",
                                                         "Old_T1")
  new_classification[i,]$PanCan_Last_Tertile = ifelse(Age_OnyCT.df[i,]$AGE > quantile(Age_OnyCT.df$AGE, prob = 0.66),
                                                         "Old_T3",
                                                         "Young_T3")
}
head(new_classification)  ## inspect

Age_OnyCT.df.Y <- cbind(Age_OnyCT.df.Y, new_classification)    ## create a single DF with all variables
head(Age_OnyCT.df.Y)    ## inspect

##
library(ggpubr)
ggviolin(data = Age_OnyCT.df.Y, x = 'TCGA_CancerCode', y = 'AGE', add = 'jitter', 
         color = c('PanCan_Last_Quartile'))



## lines 53 and 90 show NA for the second rowname - why?

## Continue classifying subtypes to get a unified data frame

## input remains constant
input <- as.data.frame(cbind(Age_OnyCT.df.Y$CASE_ID, 
                             Age_OnyCT.df.Y$AGE, 
                             Age_OnyCT.df.Y$CANCER_TYPE, 
                             Age_OnyCT.df.Y$TCGA_CancerCode))

## merge quartile columns to one. similarly, tertiles and subtype quartiles....
to_filter <- as.data.frame(cbind(Age_OnyCT.df.Y$Median_Classification,
                     Age_OnyCT.df.Y$Mean_Classification,
                     Age_OnyCT.df.Y$PanCan_First_Quartile,
                     Age_OnyCT.df.Y$PanCan_Last_Quartile,
                     Age_OnyCT.df.Y$PanCan_First_Tertile,
                     Age_OnyCT.df.Y$PanCan_Last_Tertile,
                     Age_OnyCT.df.Y$CancerAgeMedianCategory,
                     Age_OnyCT.df.Y$Subtype_Mean_Classification,
                     Age_OnyCT.df.Y$Subtype_Median_Classification,
                     Age_OnyCT.df.Y$Subtype_First_Quartile,
                     Age_OnyCT.df.Y$Subtype_Last_Quartile,
                     Age_OnyCT.df.Y$Subtype_First_Tertile,
                     Age_OnyCT.df.Y$Subtype_Last_Tertile))
colnames(to_filter) <- c("Median_Classification", "Mean_Classification","PanCan_First_Quartile",
                         "PanCan_Last_Quartile","PanCan_First_Tertile","PanCan_Last_Tertile",
                         "CancerAgeMedianCategory","Subtype_Mean_Classification",
                         "Subtype_Median_Classification","Subtype_First_Quartile",
                         "Subtype_Last_Quartile","Subtype_First_Tertile","Subtype_Last_Tertile")

## creating unified metadata
for (i in 1:nrow(to_filter)) {
  input[i,5] <- to_filter$Median_Classification[i]
  input[i,6] <- to_filter$Mean_Classification[i]
  input[i,7] <- ifelse(to_filter$PanCan_First_Quartile[i] == "Young_Q1", "Young_Q1", (ifelse(to_filter$PanCan_Last_Quartile[i] == "Old_Q3", "Old_Q3", "neither")))
  input[i,8] <- ifelse(to_filter$PanCan_First_Tertile[i] == "Young_T1", "Young_T1", (ifelse(to_filter$PanCan_Last_Tertile[i] == "Old_T3", "Old_T3", "T2")))
  input[i,9] <- to_filter$CancerAgeMedianCategory[i]
  input[i,10] <- to_filter$Subtype_Mean_Classification[i]
  input[i,11] <- to_filter$Subtype_Median_Classification[i]
  input[i,12] <- ifelse(to_filter$Subtype_First_Quartile[i] == "Young_Q1", "Young_Q1", (ifelse(to_filter$Subtype_Last_Quartile[i] == "Old_Q3", "Old_Q3", "neither")))
  input[i,13] <- ifelse(to_filter$Subtype_First_Tertile[i] == "Young_T1", "Young_T1", (ifelse(to_filter$Subtype_Last_Tertile[i] == "Old_T3", "Old_T3", "T2")))
}
colnames(input) <- c("CASE_ID","AGE","CANCER_TYPE","TCGA_CancerCode", "Median_Classification",
                     "Mean_Classification", "PanCan_Quartile", "PanCan_Tertile", "CancerAgeMedianCategory",
                     "Subtype_Mean_Classification","Subtype_Median_Classification","Subtype_Quartile", "Subtype_Tertile")

Age_OnyCT.meta <- input

#Age_OnyCT.meta
