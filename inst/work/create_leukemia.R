##  Golub data
##  Golub et al. (1999) Molecular classification of cancer: class discovery and class
##  prediction by gene expression monitoring, Science 286, 531–-537,
##  doi: 10.1126/science.286.5439.531
##
##  We take the two versions of the Leukemia data from the web site of
##  Trevor Hastie, as these data are used in the book
##  Efron, B. and Hastie, T. (2016) "Computer Age Statistical Inference".
##
##  See the data description here:
##  https://hastie.su.domains/CASI_files/DATA/leukemia.html
##
##  Another version of the complete data set (with 7029 genes) is
##      available from Bioconductor, package golubEsets, dataset Golub_Merge.
##
##  Yet another version can be found on Kaggle:
##      https://www.kaggle.com/datasets/crawford/gene-expression
##
##  Yet another version is available in the package 'plsgenomics' as
##      a dataset 'leukemia' with 3051 genes for 38 tumor mRNA samples.
##      It was taken from the package 'multtest' (archived on CRAN). See also the reference
##      S. Dudoit, J. Fridlyand and T. P. Speed (2002). Comparison of
##      discrimination methods for the classification of tumors using gene
##      expression data, Journal of the American Statistical Association 97, 77–87.
##


##  leukemia_big <- read.csv("http://hastie.su.domains/CASI_files/DATA/leukemia_big.csv")
leukemia_big <- read.csv("leukemia_big.csv")
head(leukemia_big)
leukemia_big <- t(leukemia_big)
Class <- substr(rownames(leukemia_big), 1, 3)
rownames(leukemia_big) <- NULL
leukemia_big <- as.data.frame(leukemia_big)
leukemia_big$Class <- factor(Class)
head(leukemia_big)
dim(leukemia_big)

##  leukemia_small <- read.csv("http://hastie.su.domains/CASI_files/DATA/leukemia_small.csv")
leukemia_small <- read.csv("leukemia_small.csv")
head(leukemia_small)
leukemia_small <- t(leukemia_small)
Class <- substr(rownames(leukemia_small), 1, 3)
rownames(leukemia_small) <- NULL
leukemia_small <- as.data.frame(leukemia_small)
leukemia_small$Class <- factor(Class)
head(leukemia_small)
dim(leukemia_small)

save(leukemia_big, leukemia_small, file="golub.rda")

## Retrieve the Bioconductor Golub_Merge data set and save as
##  a data frame 'golub_bioconductor' with a Class variable
library(Biobase)
library(golubEsets)
data(Golub_Merge)
x <- t(exprs(Golub_Merge))
grp <- Golub_Merge$ALL.AML

golub_bioconductor <- as.data.frame(x)
golub_bioconductor$Class <- grp
save(golub_bioconductor, file="golub_bioconductor.rda")
