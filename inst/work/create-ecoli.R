##
##  Ecoli Desease data set
##
##  UCI:    https://archive.ics.uci.edu/ml/datasets/ecoli
##
df <- read.table("ecoli.data")
rownames(df) <- df[,1]
df <- df[,-1]
Class <- factor(df[, ncol(df)])
df <- df[,-ncol(df)]
colnames(df) <- paste0("V", 1:ncol(df))
df <- cbind.data.frame(df, Class)

ecoli <- df
save(ecoli, file="ecoli.rda", version=2)
