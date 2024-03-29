if(FALSE)
{

    library(rrcov)
    library(mlrob)
    data("fruit", package="rrcov")
    obj <- LdaClassic(cultivar~., data=fruit)
    cv(obj)

    obj <- LdaClassic(Species~., data=iris)
    obj <- Linda(Species~., data=iris, l1med=TRUE)
    cv(obj)

    library(tclust)
    data(swissbank)
    head(swissbank)

    swissbank$Group <- as.factor(c(rep("genuine", 100), rep("counterfeit", 100)))

    ## obj <- LdaClassic(Group~., data=swissbank)
    obj <- Linda(Group~., data=swissbank, cov.control=CovControlMrcd(alpha=0.75))

    cv(obj)

    data(hemophilia)
    obj <- Linda(as.factor(gr)~., data=hemophilia)
    cv(obj)

    data(soil)
    soil1983 <- soil[soil$D == 0, -2]       # only 1983, remove column D (always 0)
    obj <- Linda(F~., data=soil1983, l1med=TRUE)
    cv(obj)

    data(diabetes, package="rrcov")
    obj <- Linda(group~., data=diabetes, l1med=TRUE)
    cv(obj)

    data(salmon)
    obj <- Linda(Origin~., data=salmon, l1med=TRUE)
    cv(obj)

    data(pottery)
    obj <- Linda(origin~., data=pottery, l1med=TRUE)
    cv(obj)

    data(olitos)
    obj <- Linda(grp~., data=olitos, l1med=TRUE)
    cv(obj)

    data(fish)
    # remove observation #14 containing missing value
    fish <- fish[-14,]
    # The height and width are calculated as percentages
    #   of the third length variable
    fish[,5] <- fish[,5]*fish[,4]/100
    fish[,6] <- fish[,6]*fish[,4]/100
    obj <- Linda(Species~., data=fish, l1med=TRUE)
    cv(obj)


}

if(FALSE) {

#   iris
    xx <- get_data("Iris")
    obj <- LdaFisher(xx$x, xx$grp)
    holdout(obj)
    cv(obj)
    loocv(obj)

    obj <- LdaPca(xx$x, xx$grp, k=2, preprocess="sphere")
    holdout(obj)
    cv(obj)
    loocv(obj)

## wine data
    xx <- get_data("Wine")
    obj <- LdaFisher(xx$x, xx$grp)
    holdout(obj)
    cv(obj)
    loocv(obj)

    obj <- LdaPca(xx$x, xx$grp, k=2, preprocess="sphere")
    holdout(obj)
    cv(obj)
    loocv(obj)

## diabetes data
    xx <- get_data("Diabetes")
    obj <- LdaFisher(xx$x, xx$grp)
    holdout(obj)
    cv(obj)
    loocv(obj)

    obj <- LdaPca(xx$x, xx$grp, k=2, preprocess="sphere")
    holdout(obj)
    cv(obj)
    loocv(obj)

## NIR
    xx <- get_data("NIR")
    obj <- LdaPca(xx$x, xx$grp, k=2)
    holdout(obj)
    cv(obj)

    xx <- get_data("MLL")
    obj <- LdaPca(xx$x, xx$grp, k=2)
    cv(obj)

}

cv <- function(obj, k=10){

    accur <- function(actual, predicted) {
        tab <- table(actual, predicted)
        sum(tab[row(tab)==col(tab)])/sum(tab)
    }

    is_singular <- function (mat) {
        p <- ncol(mat)
        if (!is.qr(mat))
            mat <- qr(mat)
        return(mat$rank < p)
    }

    if(missing(obj))
        return(list(aveacc=NA, acc=NA, aveaccrob=NA, accrob=NA, eaer=NA, eaerrob=NA))

    if(typeof(obj) == "S4") {
        data <- obj@X
        grp <- obj@grp
        ng <- length(levels(grp))
        method <-  if("method" %in% slotNames(class(obj))) obj@method else class(obj)
    } else {
        data <- obj$X
        grp <- obj$grp
        ng <- length(levels(grp))
        method <-  if("method" %in% names(obj)) obj$method else class(obj)
    }

    glev <- names(table(grp))
    ng <- length(glev)
    nj <- as.numeric(table(grp))
    k <- min(min(table(grp)), k)           # number of folds not less than number of obs in a class
    folds <- bfolds(grp, k)

    acc <- rep(NA, k)    # accuracies
    accrob <- rep(NA, k) # accuracies robust

    for(kk in 1:k) { # one fold

        ind <- folds[[kk]]
        dtest <- data[ind, ]
        gtest <- grp[ind]
        dtrain <- data[-ind,]
        gtrain <- grp[-ind]

        #resj <- get(dmethod)(dtrain, gtrain, method=method)
        resj <- if(is(obj, "LdaClassic")) {
                    rrcov::LdaClassic(dtrain, gtrain)
                } else if(is(obj, "Linda")){
                    rrcov::Linda(dtrain, gtrain, method=method, l1med=obj@l1med, cov.control=obj@control)
                } else if(is(obj, "QdaClassic")){
                    rrcov::QdaClassic(dtrain, gtrain)
                } else if(is(obj, "QdaCov")){
                    rrcov::QdaCov(dtrain, gtrain, method=obj@control)
                } else if(is(obj, "SosDiscRobust")){
                    rrcovHD::SosDiscRobust(dtrain, gtrain, lambda=obj@lambda)
                } else if(is(obj, "LdaFisher")){
                    LdaFisher(dtrain, gtrain)
                } else if(is(obj, "LdaPca")){
                    LdaPca(dtrain, gtrain, k=obj$k, preprocess=obj$preprocess)
                } else if(is(obj, "LdaRotatedCS")){
                    LdaRotatedCS(dtrain, gtrain, k=obj$k, preprocess=obj$preprocess)
                } else if(is(obj, "LdaPenalizedCS")){
                    LdaPenalizedCS(dtrain, gtrain, k=obj$k, preprocess=obj$preprocess, prednorm=obj$prednorm)
                } else if(is(obj, "LdaPenalizedCSV2")){
                    LdaPenalizedCSV2(dtrain, gtrain, k=obj$k, preprocess=obj$preprocess, prednorm=obj$prednorm)
                } else if(is(obj, "LdaSparse")){
                    LdaSparse(dtrain, gtrain, lambda=obj$lambda, numVars=obj$numVars, maxRSS=obj$maxRSS)
                } else if(is(obj, "LdaNSC")){
                    LdaNSC(dtrain, gtrain, threshold=obj$threshold)
                } else if(is(obj, "LdaRegularized")){
                    LdaRegularized(dtrain, gtrain, alpha=obj$alpha, delta=obj$delta, regularization=obj$regularization)
                } else {
                    stop("ERROR: unknown class")
                }

        resjpred <- predict(resj, dtest)
        gpred <- if(is(obj, "LdaFisher") ||
                    is(obj, "LdaPca") ||
                    is(obj, "LdaRotatedCS") ||
                    is(obj, "LdaPenalizedCS") ||
                    is(obj, "LdaPenalizedCSV2") ||
                    is(obj, "LdaSparse") ||
                    is(obj, "LdaNSC") ||
                    is(obj, "LdaRegularized")) resjpred$grp else resjpred@classification

        acc[kk] <- accur(gtest, gpred)

        ## Compute mahalanobis distances, using T and C from the
        ##  corresponding classifier and use them for outlier
        ##  detection (crit=qchisq...) - and then compute the
        ##  robust misclassification error
        if(typeof(resj) == "S4" && "cov" %in% slotNames(class(resj)))
        {
            md2 <- rep(NA, length(gtest))
            for(j in 1:ng){
                xcov <- if(is(obj, "QdaClassic") | is(obj, "QdaCov")) resj@cov[,,j] else resj@cov
                xinv <- if(is_singular(xcov))  solve(xcov) else pracma::pinv(xcov)

                md2[gtest==glev[j]] <- mahalanobis(dtest[gtest==glev[j],], resj@center[j,], xinv, inverted=TRUE)
            }

            crit <- qchisq(0.975, ncol(dtest))
    ##        print(md2)
    ##        print(summary(md2))
    ##        print(length(md2))
    ##        print(crit)
    ##        print(length(which(md2<crit)))

            accrob[kk] <- accur(gtest[md2 < crit], resjpred@classification[md2 < crit])
        } else {
            accrob <- NA
        }
  }

  list(aveacc=mean(acc), acc=acc, aveaccrob=mean(accrob[is.finite(accrob)]), accrob=accrob,
    eaer=1-mean(acc), eaerrob=1-mean(accrob[is.finite(accrob)]))
}

holdout <- function(obj, nsim=10, seed, test_size=0.33)
{
    accur <- function(actual, predicted) {
        tab <- table(actual, predicted)
        sum(tab[row(tab)==col(tab)])/sum(tab)
    }

    if(missing(obj))
        return(list(aveacc=NA, acc=NA, eaer=NA))

    if(typeof(obj) == "S4") {
        data <- obj@X
        grp <- obj@grp
        ng <- length(levels(grp))
        method <-  if("method" %in% slotNames(class(obj))) obj@method else class(obj)
    } else {
        data <- obj$X
        grp <- obj$grp
        ng <- length(levels(grp))
        method <-  if("method" %in% names(obj)) obj$method else class(obj)
    }
    data1 <- cbind.data.frame(data, grp=grp)

    if(!missing(seed))
        set.seed(seed)

    acctrain <- acctest <- rep(NA, nsim)    # accuracies

    for (iii in 1:nsim) {

        ## Split the data into train and test

        dtrain <- dtest <- NULL

        ## Resample
        df <- data1[sample(1:nrow(data1), nrow(data1)), ]

        for(i in levels(df$grp))
        {
            x <- df[df$grp == i, ]
            n <- nrow(x)
            ntrain <- n - round(n - (1-test_size) * n)
            ind <- sample(1:n, ntrain)
            xtrain <- x[ind,]
            xtest <- x[-ind,]
            dtrain <- rbind(dtrain, xtrain)
            dtest <- rbind(dtest, xtest)
        }

        gtrain <- dtrain[, ncol(dtrain)]
        dtrain <- dtrain[, -ncol(dtrain)]
        gtest <- dtest[, ncol(dtest)]
        dtest <- dtest[, -ncol(dtest)]

        resj <- if(is(obj, "LdaClassic")) {
                    rrcov::LdaClassic(dtrain, gtrain)
                } else if(is(obj, "Linda")){
                    rrcov::Linda(dtrain, gtrain, method=method, l1med=obj@l1med, cov.control=obj@control)
                } else if(is(obj, "QdaClassic")){
                    rrcov::QdaClassic(dtrain, gtrain)
                } else if(is(obj, "QdaCov")){
                    rrcov::QdaCov(dtrain, gtrain, method=obj@control)
                } else if(is(obj, "SosDiscRobust")){
                    rrcovHD::SosDiscRobust(dtrain, gtrain, lambda=obj@lambda)
                } else if(is(obj, "LdaFisher")){
                    LdaFisher(dtrain, gtrain)
                } else if(is(obj, "LdaPca")){
                    LdaPca(dtrain, gtrain, k=obj$k, preprocess=obj$preprocess)
                } else if(is(obj, "LdaRotatedCS")){
                    LdaRotatedCS(dtrain, gtrain, k=obj$k, preprocess=obj$preprocess)
                } else if(is(obj, "LdaPenalizedCS")){
                    LdaPenalizedCS(dtrain, gtrain, k=obj$k, preprocess=obj$preprocess, prednorm=obj$prednorm)
                } else if(is(obj, "LdaPenalizedCSV2")){
                    LdaPenalizedCSV2(dtrain, gtrain, k=obj$k, preprocess=obj$preprocess, prednorm=obj$prednorm)
                } else if(is(obj, "LdaSparse")){
                    LdaSparse(dtrain, gtrain, lambda=obj$lambda, numVars=obj$numVars, maxRSS=obj$maxRSS)
                } else if(is(obj, "LdaNSC")){
                    LdaNSC(dtrain, gtrain, threshold=obj$threshold)
                } else if(is(obj, "LdaRegularized")){
                    LdaRegularized(dtrain, gtrain, alpha=obj$alpha, delta=obj$delta, regularization=obj$regularization)
                } else {
                    stop("ERROR: unknown class")
                }

        resjpred <- predict(resj, dtest)
        gpred <- if(is(obj, "LdaFisher") ||
                    is(obj, "LdaPca") ||
                    is(obj, "LdaRotatedCS") ||
                    is(obj, "LdaPenalizedCS") ||
                    is(obj, "LdaPenalizedCSV2") ||
                    is(obj, "LdaSparse") ||
                    is(obj, "LdaNSC") ||
                    is(obj, "LdaRegularized")) resjpred$grp else resjpred@classification
        acctest[iii] <- accur(gtest, gpred)

        resjpred <- predict(resj, dtrain)
        gpred <- if(is(obj, "LdaFisher") ||
                    is(obj, "LdaPca") ||
                    is(obj, "LdaRotatedCS") ||
                    is(obj, "LdaPenalizedCS") ||
                    is(obj, "LdaPenalizedCSV2") ||
                    is(obj, "LdaSparse") ||
                    is(obj, "LdaNSC") ||
                    is(obj, "LdaRegularized")) resjpred$grp else resjpred@classification
        acctrain[iii] <- accur(gtrain, gpred)
    }

    list(acctrain=mean(acctrain), acctest=mean(acctest),
        eaertrain=round(100*(1-mean(acctrain)), 1), eaertest=round(100*(1-mean(acctest)), 1))
}

## Internal function to perform leaving-one-out cross validation by brute force -
##  recalculates the estimator n times, excluding each observation in turn.
##
##  - The discriminantion method (QDA or LDA) is selected according to the type of
##      the object.
##  - The method for computing the covariance matrices (or the common
##      cov matrix in LDA) is selected according the slot methodID.
##
loocv <- function(obj, X, grouping){

    if(missing(X))
        X <- obj$X
    if(missing(grouping))
        grp <- obj$grp
    ng <- length(levels(grp))

    ptm <- proc.time()

    n <- nrow(X)
    p <- ncol(X)
    grpnew <- rep(NA, n)

    for(i in 1:n)
    {
      ##    cat("i=",i,"\n")

        ll <- if(is(obj, "LdaFisher")) {
                LdaFisher(X[-i,], grouping=grp[-i])
            } else if(is(obj, "LdaPca")){
                LdaPca(X[-i,], grp[-i], k=obj$k, preprocess=obj$preprocess)
            } else if(is(obj, "LdaRotatedCS")){
                LdaRotatedCS(X[-i,], grouping=grp[-i], k=obj$k, preprocess=obj$preprocess)
            } else if(is(obj, "LdaPenalizedCS")){
                LdaPenalizedCS(X[-i,], grouping=grp[-i], k=obj$k, preprocess=obj$preprocess, prednorm=obj$prednorm)
            } else if(is(obj, "LdaPenalizedCSV2")){
                LdaPenalizedCSV2(X[-i,], grouping=grp[-i], k=obj$k, preprocess=obj$preprocess, prednorm=obj$prednorm)
            } else if(is(obj, "LdaSparse")){
                LdaSparse(X[-i,], grouping=grp[-i], lambda=obj$lambda, numVars=obj$numVars, maxRSS=obj$maxRSS)
            } else if(is(obj, "LdaNSC")){
                LdaNSC(X[-i,], grouping=grp[-i], threshold=obj$threshold)
            } else if(is(obj, "LdaRegularized")){
                LdaRegularized(X[-i,], grouping=grp[-i], alpha=obj$alpha, delta=obj$delta, regularization=obj$regularization)
            } else if(is(obj, "lda")){
                lda(X[-i,], grouping=grp[-i])
            } else {
                stop("ERROR: unknown class")
            }

        pp <- predict(ll, newdata=X[i,,drop=FALSE])

        grpnew[i] <- if(is(obj, "lda")) pp$class[1] else pp$grp[1]

        ##  cat("i=",i,"grp=",grp[i], "grpnew=", grpnew[i], "\n")
    }

    tab <- table(grp, grpnew)
    acc <- sum(tab[row(tab)==col(tab)])/sum(tab)

    list(tab=tab, acc=acc, eaer=1-acc)
}
