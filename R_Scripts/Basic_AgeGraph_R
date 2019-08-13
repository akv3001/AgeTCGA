setwd ("/Users/romypichardo/Downloads/cBioportal_Survival_Data")
my_data <- read.delim(file.choose()) #open file 

head(my_data)
tail(my_data)
str(my_data) ##condensed view
dim(my_data)
names(my_data)

my_data [,c("CASE_ID", "AGE","CANCER_TYPE")]   ##print column AGE AND CANCER TYPE 
my_data1 <- na.omit(my_data$AGE) ##remove NA's

Age <- my_data1  ## assign variable Age to my_data1
hist(Age, main= "Acc Age Distribution ", border = "black", col = "cadetblue2") ## plot column age as a Histogram. Main is the name of the plot. Border and color for histogram. 


