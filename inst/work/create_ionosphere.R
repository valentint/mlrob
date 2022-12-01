##
##  https://archive.ics.uci.edu/ml/datasets/ionosphere
##
##  Download the data from UCI

df <- read.csv("ionosphere.data", header=FALSE)
head(df)

grp <- df[, ncol(df)]
df <- df[, -c(1,2,ncol(df))]
head(df)
dim(df)

grp[grp=="g"] <- "good"
grp[grp=="b"] <- "bad"
grp <- factor(grp)

colnames(df) <- paste0("V", 1:ncol(df))
ionosphere <- cbind.data.frame(df, Class=grp)