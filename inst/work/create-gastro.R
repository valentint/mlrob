##  Gastrointestinal data
##
##  https://archive.ics.uci.edu/ml/datasets/Gastrointestinal+Lesions+in+Regular+Colonoscopy
##
##  http://www.depeca.uah.es/colonoscopy_dataset/
##
##  The file contains the following info:
##  1st row: name of the lesion
##  2nd row: class of lesion
##  3rd row: type of light used (one column is with White Light (WL) and the other with Narrow Band Imaging (NBI).
##      This makes that, instead of only 76 columns we have 152)
##  4rt row until the last one: features

df <- read.csv("gastro.txt")

## First row: class

## Second row: light
df1 <- df[, df[2,] == 1]
df2 <- df[, df[2,] == 2]
colnames(df1) <- NULL
colnames(df2) <- NULL

grp <- df1[1,]
rownames(grp) <- NULL

df1 <- df1[-c(1:2), ]
df2 <- df2[-c(1:2), ]

dim(df1)

df <- cbind.data.frame(t(df1), t(df2))
colnames(df) <- paste0("V", 1:ncol(df))
df <- cbind.data.frame(class=t(grp), df)

gastro <- df
save(gastro, file="gastro.rda")
