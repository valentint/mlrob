if(FALSE) {
    library(pracma)
    library(mlrob)

    xx <- get_data("Iris")
    ldar <- LdaPenalizedCSV2(xx$x, xx$grp, preprocess="sphere", prednorm=TRUE)

    (tab <- table(xx$grp, predict(ldar, newdata=xx$x)$grp))
    round(100*(1-sum(diag(tab))/sum(tab)),1)

    round(100*(1-cv(ldar)$aveacc),1)
    round(100*loocv(ldar)$eaer,1)
    holdout(ldar)
}

##  'prednorm\ is an argument used for prediction: whether to normalize the scores
##
LdaPenalizedCSV2 <- function(x, grouping, prior=proportions, k=ncol(x),
    preprocess=c("center", "sphere", "standardize"), tol=1e-6, prednorm=FALSE, trace=FALSE) {

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
         else if(preprocess == "sphere") {
            x1 <- scale(x, center=TRUE, scale=FALSE)
            vx <- apply(x, 2, function(x) sd(x) / sqrt(1/(n-1)))
            id <- which(vx > 0)
            z1 <- scale(x1[, id], center=FALSE, scale=vx[id])
            attr(z1, "scaled:center") <- attr(x1, "scaled:center")
            attr(z1, "scaled:scale") <- vx
            z1
         } else if(preprocess == "standardize") {
            x1 <- scale(x, center=TRUE, scale=FALSE)
            vx <- apply(x, 2, sd)
            id <- which(vx > 0)
            z1 <- scale(x1[, id], center=FALSE, scale=vx[id])
            attr(z1, "scaled:center") <- attr(x1, "scaled:center")
            attr(z1, "scaled:scale") <- vx
            z1
         } else
            stop("Wrong preprocessing selected!")

    scaled_center <- if(!is.null(attr(z, "scaled:center"))) attr(z, "scaled:center") else rep(0, p)
    scaled_scale <- if(!is.null(attr(z, "scaled:scale"))) attr(z, "scaled:scale") else rep(1, p)

    z <- as.data.frame(z)       # to delete the attributes
    z <- as.matrix(z)
    p <- ncol(z)
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
    F0 <- F <- W$u
    per <- W$d

    if(trace)
        cat("\nPercent variance explained:", 100*sum(per[1:r]^2)/sum(per^2), "\n")

    per <- per[1:r]
    F <- F[,1:r]

    ## Keep the PCA parameters
    loadings <- W$v[, 1:r]
    sdev <- W$d[1:r]/sqrt(max(1, n - 1))

##    F1 <- z %*% loadings %*% diag(1/per)
##    F2 <- z %*% loadings %*% diag(1/(sdev*sqrt(max(1, n - 1))))

##    print(head(F))
##    print(head(F1))
##    print(head(F2))

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

    ##V <- eigen(t(F) %*% (H + zz) %*% F)
    ##Q <- V$vectors
    ##d <- V$values

    rg <- r                 # actual dimensions
    if(n > p) {
    	rg <- sum(d > tol)
        rg <- min(rg, r)
    }

    if(trace)
        cat("\nNumber of actual dimensions used: ", rg, "\n")

    F0 <- F0[, 1:rg]

    ##Q <- Q[, 1:rg]
    ##F <- F %*% Q

    F <- F[, 1:rg]              # truncated scores
    Cx <- C <- W %*% F                # centroids of the scores

    W <- svd(t(F) %*% F0)
    Q <- W$v %*% t(W$u)

    rotation <- Q

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
                meanj=meanj, cv=cv, meanov=meanov, nobs=n, centroids=Cx,       #    scores=fs, fdiscr=fdiscr,
                sdev=sdev, loadings=loadings, rotation=rotation,
                SB=SB, SB_final=SBN,
                mc=mc, mcrate=rate, grppred=grppred,
                method="LDA on penalized component scores",
                X=x, grp=g, k=r, rg=rg, prednorm=prednorm,
                preprocess=preprocess, scaled_center=scaled_center, scaled_scale=scaled_scale)
    class(res) <- "LdaPenalizedCSV2"
    res
}

print.LdaPenalizedCSV2 <- function(x,...){
  cat("--------------------------------------")
  cat("\nResults from Fishers LDA on penalized components")
  cat("\n- Initial between-groups scatter...:", sum(diag(x$SB)))
  cat("\n- Resulting between-grops scatter.: ", sum(diag(x$SB_final)), "\n")

  ##    cat("\n- Loadings matrix: \n")
  print(x$loadings)
  cat("--------------------------------------\n")
}

predict.LdaPenalizedCSV2 <- function(object, newdata){
    ct <- FALSE
    if(missing(newdata)) {
        newdata <- object$X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)
    n <- nrow(x)
    p <- ncol(x)

    ## Normalize and predict
    id <- which(object$scaled_scale > 0)
    z <- scale(x[, id, drop=FALSE], center=object$scaled_center[id], scale=object$scaled_scale[id])

    if(length(object$meanj) > 0 & ncol(z) != nrow(object$meanj) | ncol(z) != nrow(object$loadings))
        stop("wrong number of variables")

    ng <- ncol(object$meanj) # number of groups

    ## Get the initial scores using 'loadings' and 'sdev'
    if(!object$prednorm) {
        F1 <- z %*% object$loadings %*% diag(1/(object$sdev*sqrt(max(1, object$nobs - 1))))
        F1 <- F1 %*% object$rotation
    } else {
        F1 <- z %*% object$loadings
        F1 <- F1 %*% object$rotation
        dVV <- diag(t(F1) %*% F1)
        F1 <- F1 %*% diag(1/sqrt(dVV))              # = ZA
    }

    # Restore and rotate the centroids
    C1 <- t(object$meanj) %*% object$loadings %*% diag(1/(object$sdev*sqrt(max(1, object$nobs - 1))))
    C1 <- C1 %*% object$rotation

    ## Discriminant scores
    fs <- matrix(NA, nrow=n, ncol=ng)
    for(j in 1:ng) {
        xc <- scale(F1, object$centroids[j,], scale=FALSE)
        fs[,j] <- sqrt(apply(xc^2, 1, sum))
    }

    ## prediction:
    grp <- apply(fs, 1, which.min)
    grpnam <- levels(object$grp)[grp]
    grp <- factor(grpnam, levels=levels(object$grp))

    ret <- list(grp=grp)
    if(ct) {
        ret$ct <- table(object$grp, grp)
        ret$AER <- 1 - sum(diag(ret$ct))/sum(ret$ct)
    }

    ret
}
