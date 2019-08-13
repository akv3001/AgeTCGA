setwd("/Users/romypichardo/Downloads/cBioportal_test") ## set directory

file = list.files(pattern="\\.txt$") ## assign all files to variable file 

list_clinical = list() ## create emapty list with list()

list_clinical[["Cancer Type"]] = input_read[,c("AGE")]
