######
##  VT::15.02.2023
##
##
##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/R/mlrob", load_code=roxygen2:::load_installed)
##
#' Balanced folds for cross-validation for classification
#'
#' @description Generates balanced folds for cross-validation for classification
#'
#' @param y a factor with class labels
#' @param nfolds number of folds, defaults to 10 but should not be
#'  larger than the number of objects in the smallest class
#'
#' @return  a list of \code{nfolds} vectors of length approx \code{length(y)/nfolds} each
#'
#' @examples
#'
#'  folds <- bfolds(iris$Species)
#'  length(folds)
#'  folds[[1]]
#'
#' @export
#' @author Valentin Todorov \email{valentin@@todorov.at}
bfolds <- function (y, nfolds=min(min(table(y)), 10))
{

    permute.rows <- function (x){
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
    }

    totals <- table(y)
    fmax <- max(totals)
    nfolds <- min(nfolds, fmax)
    folds <- as.list(seq(nfolds))
    yids <- split(seq(y), y)
    bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
    for (i in seq(totals)) {
        bigmat[seq(totals[i]), i] <- sample(yids[[i]])
    }
    smallmat <- matrix(bigmat, nrow = nfolds)
    smallmat <- permute.rows(t(smallmat))
    res <- vector("list", nfolds)
    for (j in 1:nfolds) {
        jj <- !is.na(smallmat[, j])
        res[[j]] <- smallmat[jj, j]
    }
    return(res)
}
