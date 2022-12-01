##
##  BUPA
##
##  https://archive.ics.uci.edu/ml/datasets/liver+disorders
##
df <- read.csv("bupa.data", header=FALSE)
head(df)

colnames(df) <- c("mcv", "alkphos", "sgpt", "sgot", "gammagt", "drinks", "selector")
bupa_original <- df
bupa <- df[, 1:5]
bupa$Class <- factor(ifelse(df[,6] < 5, 1, 2))

##  Remove four duplicated rows
##
##  Thanks to Leon for mentioning that there are duplicates in this data set.
##  --UCI ML Librarian

##  row 84 and 86:   94,58,21,18,26,2.0,2
##  row 141 and 318:   92,80,10,26,20,6.0,1
##  row 143 and 150:   91,63,25,26,15,6.0,1
##  row 170 and 176:   97,71,29,22,52,8.0,1

bupa <- bupa[-c(86, 318, 150, 176),]

save(bupa, bupa_original, file="bupa.rda")
