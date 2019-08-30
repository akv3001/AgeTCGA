setwd("/Users/romypichardo/Desktop/BigData_Project/AgeTCGA/RNASeq_Differential_Exp/Samples")
path = ("RNASeq_Differential_Exp/Samples")
list.files("RNASeq_Differential_Exp/Samples")


data <- data.matrix(read.table("Counts_OVEREXP_XBP1s_input.txt", row.names = 1, header = T))

# Review dataset
dim(data)    ##size of dataset
head(data)   ##First 6 rows
tail(data)   ##Last 6 rows

## Overview plot
hist(data, col = "blue", main="Sample44")

install.packages("rgl")
library(rgl)
plot3d(pc$scores[,1:3], col=iris$Species)





##data1 <- na.omit(data) ## deletion of missing
##data1 <- scale(data1) ##standardize variables


##log 2 plot results
data2 = log2(data)
# View data after log-transformation
hist(data2, col = "blue", main="Sample44 (log2)")

# Separate the two conditions into two smaller data frames
X = data2[,1:3]
C = data2[,4:6]


# Calculate the means for each sample group 
X.mean = apply(X, 1, mean)
C.mean = apply(C, 1, mean)

########
# maximum of all the means
limit = max(X.mean, C.mean)


# Scatter plot mean
plot(C.mean ~ X.mean, xlab = "X", ylab = "C",
     main = "Sample44- ScatterPlot", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "RED")

# Compute fold-change (biological significance)
# Difference between the means of the conditions
fold = X.mean - C.mean

#### Histogram of the fold differences
hist(fold, col = "blue")

