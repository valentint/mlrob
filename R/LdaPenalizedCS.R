#' @importFrom e1071 matchClasses
require(e1071)

if(FALSE) {

## iris data
    grp <- iris$Species
    x <- iris[, 1:4]

    ldar <- LdaPenalizedCS(x, grp, preprocess="center")

    (tab <- table(grp, predict(ldar, newdata=x)$grp)); rrcov:::.AER(tab)
    1-cv(ldar)$aveacc
    loo(ldar)$eaer
    holdout(ldar)

## diabetes data
    data(diabetes, package="rrcov")
    x <- diabetes[, c("glucose", "insulin", "sspg")]
    grp <- diabetes$group

    ldar <- LdaRotatedCS(x, grp)

    (tab <- table(grp, predict(ldar, newdata=x)$grpnam)); rrcov:::.AER(tab)
    1-cv(ldar)$aveacc
    loo(ldar)$eaer
    holdout(ldar)

## wine data
    library(MBCbook)
    data(wine27)

    x <- wine27[, 1:(ncol(wine27)-2)]
    grp <- wine27$Type

    ldar <- LdaRotatedCS(x, grp)

    (tab <- table(grp, predict(ldar, newdata=x)$grpnam)); rrcov:::.AER(tab)
    1-cv(ldar)$aveacc
    loo(ldar)$eaer
    holdout(ldar)
}

LdaPenalizedCS <- function(x, grouping, prior=proportions, k=ncol(x), preprocess=c("center", "sphere", "standardize"), tol=1e-6, trace=FALSE) {

    preprocess <- match.arg(preprocess)

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
    } else if(length(grouping) > 1 && length(grouping) < n) {
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

    if(trace) {
        cat("\nRank of X: ", r, robustbase::rankMM(x), "\n")
        cat("\nNumber of groups: ", ng, "\n")
    }

    gx <- as.numeric(g)     # this will be used throughout the procedure

    z <- if(preprocess == "center") scale(x, center=TRUE, scale=FALSE)
         else if(preprocess == "sphere") scale(x, center=TRUE, scale=apply(x, 2, function(x) sd(x) / sqrt(1/(n-1))))
         else if(preprocess == "standardize") scale(x, center=TRUE, scale=TRUE)
         else
            stop("Wrong preprocessing selected!")

    scaled_center <- if(!is.null(attr(z, "scaled:center"))) attr(z, "scaled:center") else rep(0, p)
    scaled_scale <- if(!is.null(attr(z, "scaled:scale"))) attr(z, "scaled:scale") else rep(1, p)

    z <- as.data.frame(z)       # to delete the attributes
    z <- as.matrix(z)
    zz <- z %*% t(z)

    r <- if(n > p) p else n - 1     # rank of X
    rr <- robustbase::rankMM(x)
    r <- min(r, rr)
    r <- min(r, k)                  # number of initial components

    ## priors, means and covariances of the groups
    meanj <- matrix(NA, nrow=p, ncol=ng)
    cv <- list()
    for(j in 1:ng){
        meanj[,j] <- apply(z[g == lev1[j], ], 2, mean)
        cv[[j]] <- cov(z[g == lev1[j], ])
    }
    meanov <- meanj %*% prior

    ## Initial component scores
    W <- svd(z)
    F <- W$u
    per <- W$d

    if(trace)
        cat("\nPercent variance explained:", 100*sum(per[1:r]^2)/sum(per^2), "\n")

    per <- per[1:r]
    s <- r + ng
    F <- F[,1:r]

    ## Initial indicator matrix G
    G <- matrix(0, nrow=n, ncol=ng)
    G[(1:n) + n*(unclass(gx)-1)] <- 1               #if g is factor, the expression is G[(1:n) + n*(unclass(gx)-1)] <- 1

    ## Centroids and related stuff
    W <- diag(1/counts) %*% t(G)
    C <- W %*% F
    H <- G %*% W

    SB <- t(F) %*% H %*% F
    if(trace)
        cat("\nInitial between-grops scatter: ", sum(diag(SB)), "\n")

    ## Better scores F
    V <- eigen(H + zz)
    F <- V$vectors
    d <- V$values

    rg <- r                 # actual dimensions
    if(n > p) {
    	rg <- sum(d > tol)
    }

    if(trace)
        cat("\nNumber of actual dimensions used: ", rg, "\n")

    Fx <- F <- F[, 1:rg]              # truncated scores
    Cx <- C <- W %*% F                # centroids of the scores


    ##  Prediction using F and C
    ## Fisher scores
    fs <- matrix(NA, nrow=n, ncol=ng)
    for(j in 1:ng) {
        xc <- scale(F, C[j,], scale=FALSE)
        fs[,j] <- sqrt(apply(xc^2, 1, sum))
    }
    grppred <- apply(fs, 1, which.min)
    mc <- table(grppred, gx)
    rate <- 1 - sum(mc[row(mc) == col(mc)])/sum(mc)

    if(trace) {
        cat("\nConfusion matrix: \n")
        print(mc)
        cat("\nMisclassified observations: ", sum(mc) - sum(diag(mc)), "\n")
        cat("\nPercent misclassified: ", round(100*rate, 2), "\n")
    }

    ## Centroids and related stuff
    counts1 <- colSums(G)
    W <- diag(1/counts1) %*% t(G)
    C <- W %*% F
    H <- G %*% W
    SBN <- t(F) %*% H %*% F
    if(trace)
        cat("\nResulting between-grops scatter: ", sum(diag(SBN)), "\n")

##  ==============================================================

    res <- list(call=xcall, counts=counts,
                meanj=meanj, cv=cv, meanov=meanov, nobs=n,       #    scores=fs, fdiscr=fdiscr,
#                sdev=sdev, loadings=loadings, rotation=rotation,
                scores=Fx, centroids=Cx,
                SB=SB, SB_final=SBN,
                mc=mc, mcrate=rate, grppred=grppred,
                method="LDA on penalized component scores",
                X=x, grp=g, k=r, rg=rg,
                preprocess=preprocess, scaled_center=scaled_center, scaled_scale=scaled_scale)
    class(res) <- "LdaPenalizedCS"
    res
}

print.LdaPenalizedCS <- function(x,...){
  cat("--------------------------------------")
  cat("\nResults from Fishers LDA on penalized components")
  cat("\n- Initial between-groups scatter...:", sum(diag(x$SB)))
  cat("\n- Resulting between-grops scatter.: ", sum(diag(x$SB_final)), "\n")

  ##    cat("\n- Loadings matrix: \n")
  print(x$loadings)
  cat("--------------------------------------\n")
}

predict.LdaPenalizedCS <- function(object, newdata){
    ct <- FALSE
    if(missing(newdata)) {
        newdata <- object$X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)
    n <- nrow(x)
    p <- ncol(x)
    z <- scale(x, center=object$scaled_center, scale=object$scaled_scale)
    z1 <- scale(object$X, center=object$scaled_center, scale=object$scaled_scale)

    if(length(object$meanj) > 0 & ncol(x) != nrow(object$meanj))
        stop("wrong number of variables")

    ng <- ncol(object$meanj) # number of groups

    Fx <- object$scores
    Cx <- object$centroids

    ## According to the notes on page 16
    ##
    ##  Fx and Cx are the penalized scores and their centroids
    ##  obtained on the training sample
    D <- diag(t(Fx) %*% z1 %*% t(z1) %*% Fx)
    Fnew <- z %*% t(z1) %*% Fx %*% diag(1/D)
    Fnew <- Fnew %*% diag(1/sqrt(diag(t(Fnew) %*% Fnew)))       # normalize Fnew to unit length

    ## Discriminant scores
    fs <- matrix(NA, nrow=n, ncol=ng)
    for(j in 1:ng) {
        xc <- scale(Fnew, Cx[j,], scale=FALSE)
        fs[,j] <- sqrt(apply(xc^2, 1, sum))
    }

    ## prediction:
    grp <- apply(fs, 1, which.min)
    grpnam <- levels(object$grp)[grp]
    grp <- factor(grp, levels=levels(object$grp))

    ret <- list(grpnam=grpnam, grp=grp)
    if(ct) {
        ret$ct <- rrcov:::mtxconfusion(object$grp, grp)
        ret$AER <- 1 - sum(ret$ct[row(ret$ct) == col(ret$ct)])/sum(ret$ct)
    }

    ret
}
