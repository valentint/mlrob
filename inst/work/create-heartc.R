##
##  Heart Desease data set
##
##  UCI:    https://archive.ics.uci.edu/ml/datasets/heart+disease
##
df <- read.csv("processed.cleveland.data", header=FALSE, na.string="?")
df[is.na(df[,12]),]
df[is.na(df[,13]),]
table(df[,14])
heartc <- df

save(heartc, file="heartc.rda", version=2)

## 1. There are 6 missing values in V12 and V13
## 2. The grouping variable V14 has 5 values: recode it as 0 and 1 (if > 0)
##
##    0   1   2   3   4
##  160  54  35  35  13
