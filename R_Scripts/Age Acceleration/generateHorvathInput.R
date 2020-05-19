### This file accepts methylation beta values
### and outputs them in a format compatible with
### the Horvath clock web tool


load("~/Documents/Elemento/AgeTCGA-master/DATA/cBioportal_Survival_Data/age_metadata.RData")


PATH <- "/Users/Yajas/Documents/Elemento/AgeTCGA-master/DATA/HG 19 Meth/Beta/"
files <- list.files(PATH)
normals <- files[endsWith(files, "Normal.txt")]
files <- files[!endsWith(files, "Normal.txt")]

horvath_dir <- "/Users/Yajas/Documents/Elemento/AgeTCGA-master/DATA/HG 19 Meth/Horvath/Input/"


for (i in 1:length(files)) {
  cancer <- gsub(".*[_]([^.]+)[.].*", "\\1", files[i])
  dat0 <- read.table(paste0(PATH, files[i]))
  
  if (sum(duplicated(substr(gsub('\\.', '-', colnames(dat0)),1,15))) > 0){
    dat0 <- dat0[,!duplicated(substr(gsub('\\.', '-', colnames(dat0)),1,15))]
  }
  colnames(dat0) <- substr(gsub('\\.', '-', colnames(dat0)),1,15)
  
  common_cases <- intersect(Age_OnyCT.meta$CASE_ID, colnames(dat0))
  dat0 <- dat0[, common_cases]
  age_meta <- Age_OnyCT.meta[Age_OnyCT.meta$CASE_ID %in% common_cases,]
  
  dat0 <- dat0[,order(colnames(dat0))]
  age_meta <- age_meta[order(age_meta$CASE_ID), ]
  
  setDT(dat0, keep.rownames = TRUE)[]
  colnames(dat0)[1] <- "probeID"
  age_meta <- data.frame(Sample = age_meta$CASE_ID,
                         Age = age_meta$AGE)
  
  if(!all(age_meta$Sample == colnames(dat0)[2:ncol(dat0)])) break()
  
  write.csv(x = dat0, file = paste0(horvath_dir,cancer,"_horvath_beta.csv"), row.names = FALSE)
  write.csv(x = age_meta, file = paste0(horvath_dir,cancer,"_horvath_metadata.csv"), row.names = FALSE)
}



#### Only for normal 
cancer <- gsub(".*[_]([^.]+)[.].*", "\\1", normals)
dat0 <- read.table(paste0(PATH, normals))

filtered_age <- Age_OnyCT.meta
filtered_age$CASE_ID <- substr(gsub('\\.', '-', filtered_age$CASE_ID), 1, 12)
filtered_age <- filtered_age[!duplicated(filtered_age$CASE_ID), ]

if (sum(duplicated(substr(gsub('\\.', '-', colnames(dat0)), 1, 12))) > 0) {
  dat0 <- dat0[, !duplicated(substr(gsub('\\.', '-', colnames(dat0)), 1, 12))]
}
colnames(dat0) <- substr(gsub('\\.', '-', colnames(dat0)), 1, 12)

common_cases <- intersect(filtered_age$CASE_ID, colnames(dat0))
dat0 <- dat0[, common_cases]
age_meta <- filtered_age[filtered_age$CASE_ID %in% common_cases, ]

dat0 <- dat0[, order(colnames(dat0))]
age_meta <- age_meta[order(age_meta$CASE_ID),]

setDT(dat0, keep.rownames = TRUE)[]
colnames(dat0)[1] <- "probeID"
age_meta <- data.frame(Sample = age_meta$CASE_ID,
                       Age = age_meta$AGE)

all(age_meta$Sample == colnames(dat0)[2:ncol(dat0)])

write.csv(
  x = dat0,
  file = paste0(horvath_dir, cancer, "_horvath_beta.csv"),
  row.names = FALSE
)
write.csv(
  x = age_meta,
  file = paste0(horvath_dir, cancer, "_horvath_metadata.csv"),
  row.names = FALSE
)
