setwd("/Users/romypichardo/Desktop/BigData_Project/AgeTCGA/RNASeq_Differential_Exp/Samples")
path = ("RNASeq_Differential_Exp/Samples")


############## analysis #########
data <- data.matrix(read.table("Counts_OVEREXP_XBP1s_input.txt", row.names = 1, header = T))

# Review dataset
dim(data)    ##size of dataset
head(data)   ##First 6 rows
tail(data)   ##Last 6 rows

pdf("Sample44_Overview.pdf")
## Overview plot
hist(data, col = "blue", main="Sample44 - Overview")
dev.off()

##log 2 plot results
data2 = log2(data)

pdf("Sample44_Log_Transformation.pdf")
# View data after log-transformation
hist(data2, col = "blue", main="Sample44 Transformation(log2)")
dev.off()

# Separate the two conditions i- create two new df
X = data2[,1:3]
C = data2[,4:6]


# Calculate the means for each sample group 
X.mean = apply(X, 1, mean)
C.mean = apply(C, 1, mean)

########
# maximum of all the means
limit = max(X.mean, C.mean)


# Scatter plot mean
pdf("Sample44_ScatterPlot_Mean.pdf")
plot(C.mean ~ X.mean, xlab = "X", ylab = "C",
     main = "Sample44- ScatterPlot Mean", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "red")
dev.off()

# Compute fold-change (biological significance)
# Difference between the means of the conditions

fold = X.mean - C.mean
pdf("Sample44_Fold_Differences.pdf")
#### Histogram of the fold differences
hist(fold, col = "blue", main ="Sample44 Fold Differences" )
dev.off()



