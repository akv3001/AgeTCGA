setwd ("/Users/romypichardo/Downloads/cBioportal_Survival_Data") ## Set 

my_data <- read.delim(file.choose())
y1 <- na.omit(my_data$AGE)

my_data <- read.delim(file.choose()) #open file 
y2 <- na.omit(my_data$AGE)

my_data <- read.delim(file.choose()) #open file 
y3 <- na.omit(my_data$AGE)

my_data <- read.delim(file.choose()) #open file 
y4 <- na.omit(my_data$AGE)

##Basic BoxPlot with color selection (multi-color) strong
boxplot(y1, y2, y3, y4, 
        col=topo.colors(4), ylab="Age", 
        xlab = "Cancer Types", 
        main ="Cancer Age Distribution", 
        names=c("Thymic", "Thyroid", "Uterine Carcinosarcoma", "Glioma"))


## BoxPlot with color selection soft color strong 
boxplot(y1, y2, y3, y4,
        col=cm.colors(5), ylab="Age",
        xlab = "Cancer Types", main ="Cancer Age Distribution", 
        names=c("Non-Hodgkin Lymphoma", "Esophagogastric", "Glioblastoma Multiforme", "Uveal Melanoma"))

##BoxPlot - no color selection plain
boxplot(y1, y2, y3, y4, 
        ylab="Age", xlab = "Cancer Types", 
        main ="Cancer Age Distribution", 
        names=c("Non-Hodgkin Lymphoma", "Esophagogastric", "Glioblastoma Multiforme", "Uveal Melanoma"))
