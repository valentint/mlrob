if(FALSE) {

    ## iris data
    grp <- iris$Species
    x <- iris[, 1:4]

    (ll <- LdaNSC(x, grp))
    (pred <- predict(ll))
    1-sum(diag(pred$ct))/sum(pred$ct)

    (ll <- LdaNSC(x, grp, threshold=1))
    (pred <- predict(ll))
    1-sum(diag(pred$ct))/sum(pred$ct)



    (tt <- table(grp, predict(ll, newdata=x)$grp))             # resubstitution
    (resub <- round(100*rrcov:::.AER(tt), 1))
    (ho <- holdout(ll))                                        # Hold out validation
    (cverr <- round(100*(1-cv(ll)$aveacc), 1))                 # Cross-validation
    (looerr <- round(100*loocv(ll)$eaer, 1))                   # Leave one out cross validation
}

LdaNSC <- function(x, grouping, prior=proportions, threshold=NULL, n.threshold=30){

    if(is.null(dim(x)))
        stop("x is not a matrix")

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    if(length(grouping) == 1) {
        # this is the number of groups and the groups are of equal size
        ng = grouping
        ni = n/ng
        if(ng*ni < n)
            stop("nrow(x) is not divisible by the number of groups")
        grouping <- rep(0,0)
        for(i in 1:ng)
            grouping <- c(grouping, rep(i,ni))
    }else if(length(grouping) > 1 && length(grouping) < n) {
        # grouping contains a vector with the group sizes
        ng <- length(grouping)
        if(sum(grouping) != n)
            stop("nrow(x) is not equal to n1+n2+...+nn")

        gx <- rep(0,0)
        for(i in 1:ng)
            gx <- c(gx, rep(i,grouping[i]))
        grouping <- gx
    }

    if(n != length(grouping))
        stop("nrow(x) and length(grouping) are different")

    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))

    if(!missing(prior)) {
        if(any(prior < 0) || round(sum(prior), 5) != 1)
            stop("invalid prior")
        if(length(prior) != nlevels(g))
            stop("prior is of incorrect length")
        prior <- prior[counts > 0]
    }
    if(any(counts == 0)) {
        warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
        lev1 <- lev[counts > 0]
        g <- factor(g, levels=lev1)
        counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    names(g) <- NULL
    names(prior) <- levels(g)

    ## priors, means and covariances of the groups
    meanj <- matrix(NA, nrow=p, ncol=ng)
    cv <- list()
    for(j in 1:ng){
        meanj[,j] <- apply(x[g == lev1[j], ], 2, mean)
        cv[[j]] <- cov(x[g == lev1[j], ])
    }
    meanov <- meanj %*% prior

    mydata <- list(x=t(x), y=g)
    pp <- pamr.train(mydata, threshold=threshold, n.threshold=n.threshold)

    ## prediction:
    errmin <- which.min(pp$errors)
    threshold <- pp$threshold[errmin]
    numVars <- pp$nonzero[errmin]
    grppred <- pamr.predict(pp, mydata$x, threshold=threshold)

    ## misclassification rate:
    mc <- table(g, grppred)
    ##  mc <- mc[, matchClasses(mc, method = "exact", verbose=FALSE)]
    rate <- 1 - sum(diag(mc)) / sum(mc)

    res <- list(call=xcall, prior=prior, counts=counts,
                meanj=meanj, cv=cv, meanov=meanov,
                mc=mc, mcrate=rate, grppred=grppred,
                fit=pp, threshold=threshold, numVars=numVars,
                method="Nearest Shrunken Centroids",
                X=x, grp=g)
    class(res) <- "LdaNSC"
    res
}

print.LdaNSC <- function(x,...){
  cat("--------------------------------------")
  cat("\nResults from Nearest Shrunken Centroids")
    print(x$fit)
  cat("--------------------------------------\n")
}

predict.LdaNSC <- function(object, newdata){
    ct <- FALSE
    if(missing(newdata)) {
        newdata <- object$X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)

    if(length(object$meanj) > 0 & ncol(x) != nrow(object$meanj))
        stop("wrong number of variables")
    ng <- ncol(object$meanj) # number of groups


    ## prediction:
    grp <- pamr.predict(object$fit, t(x), threshold=object$threshold)
    grpnam <- as.character(grp)

    ret <- list(grpnam=grpnam, grp=grp)

    if(ct)
        ret$ct <- rrcov::mtxconfusion(object$grp, ret$grp)

    ret
}
