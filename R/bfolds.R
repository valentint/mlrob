## Balanced folds for cross validation
##  from rda
##
bfolds <- function (y, nfolds = min(min(table(y)), 10))
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
