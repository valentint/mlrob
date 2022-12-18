alltr <- readRDS("leukemia6train.rds")
allts <- readRDS("leukemia6test.rds")

x <- alltr[, -ncol(alltr)]
colnames(x) <- paste0("V", 1:ncol(x))

grp <- factor(alltr[, ncol(alltr)])
dim(x)
table(grp)

ALL_train <- cbind.data.frame(x, Class=grp)

x <- allts[, -ncol(allts)]
colnames(x) <- paste0("V", 1:ncol(x))
grp <- factor(allts[, ncol(allts)])
dim(x)
table(grp)

ALL_test <- cbind.data.frame(x, Class=grp)

save(ALL_train, ALL_test, file="Yeoh2002.rda", version=2)
