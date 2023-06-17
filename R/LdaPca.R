if(FALSE) {

    library(mlrob)

    xx <- get_data("Gastro")
    ldar <- LdaPca(xx$x, xx$grp)

    (tab <- table(xx$grp, predict(ldar, newdata=xx$x)$grp))
    round(100*(1-sum(diag(tab))/sum(tab)),1)
    round(100*(1-cv(ldar)$aveacc),1)
    round(100*loocv(ldar)$eaer,1)
    holdout(ldar)

}

LdaPca <- function(x, grouping, k=ncol(x), prior=proportions,
    preprocess=c("none", "center", "sphere", "standardize"), tol=1e-7, trace=FALSE){

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


    z <- if(preprocess == "none") x
         else if(preprocess == "center") scale(x, center=TRUE, scale=FALSE)
         else if(preprocess == "sphere") scale(x, center=TRUE, scale=apply(x, 2, function(x) sd(x) / sqrt(1/(n-1))))
         else if(preprocess == "standardize") scale(x, center=TRUE, scale=TRUE)
         else
            stop("Wrong preprocessing selected!")

    scaled_center <- if(!is.null(attr(z, "scaled:center"))) attr(z, "scaled:center") else rep(0, p)
    scaled_scale <- if(!is.null(attr(z, "scaled:scale"))) attr(z, "scaled:scale") else rep(1, p)

    z <- as.data.frame(z)       # to delete the attributes

    r <- if(n > p) p else n - 1     # rank of X
    rr <- robustbase::rankMM(x)
    r <- min(r, rr)
    r <- min(r, k)                  # number of initial components

    ##        cat("\nOld rank: ", r)

    pca <- prcomp(z, rank.=r)

    ##  Check for constant PCA scores
    vv <- rep(0, r)
    for(j in 1:ng){
        vvx <- apply(pca$x[g == lev1[j], ], 2, var)
        vv <- vv + vvx * proportions[j]
    }

    if(length(id <- which(vv <= tol)) > 0) {
        if((r1 <- min(id)) > 1) {
            r <- min(r, r1-1)
            cat(" New rank: ", r, "\n")

            ##  pca <- prcomp(z, rank.=r)
            pca$rotation <- pca$rotation[,1:r]
            pca$x <- pca$x[,1:r]
        }
    }

    lda <- lda(pca$x, g)

    pred <- predict(lda)
    grppred <- pred$class

    ## misclassification rate and discriminant scores
    mc <- table(g, grppred)
    rate <- 1 - sum(diag(mc)) / sum(mc)
    fdiscr <- pred$x               # discriminant scores

    res <- list(pcaobj=pca, ldaobj=lda,
                fdiscr=fdiscr,
                mc=mc, mcrate=rate, grppred=grppred,
                X=x, grp=g, k=r,
                preprocess=preprocess, scaled_center=scaled_center, scaled_scale=scaled_scale)

    class(res) <- "LdaPca"

  res
}

print.LdaPca <- function(x,...){
  cat("--------------------------------------\n")
    print(x$ldaobj)

  cat("--------------------------------------\n")
}

plot.LdaPca <- function(x, ...) {
    proj <- data.frame(x$pcaobj$x[,1:2])
    proj$grp <- as.factor(x$grp)
    proj$grppred <- as.factor(x$grppred)
    firstscores <- NULL
    secondscores <- NULL
    colnames(proj) <- c("firstscores", "secondscores","grp", "grppred")
    gg <- ggplot(proj, aes_string("firstscores", "secondscores", colour="grp", shape="grppred"))
    gg <- gg + geom_point()
    gg <- gg + xlab("First PC") + ylab("Second PC")
    gg <- gg + labs(shape="Predicted", color="Original")
    print(gg)
}

predict.LdaPca <- function(object, newdata, ...){
    ct <- FALSE
    if(missing(newdata)) {
        newdata <- object$X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)
    z <- scale(x, center=object$scaled_center, scale=object$scaled_scale)

    ng <- nrow(object$ldaobj$means)      # number of groups in 'lda' from MASS

    pr_pca <- predict(object$pcaobj, newdata=z)
    pr_lda <- predict(object$ldaobj, newdata=pr_pca)

    ret <- list(grpnam=pr_lda$class, grp=as.numeric(pr_lda$class))

    if(ct)
        ret$ct <- rrcov::mtxconfusion(object$grp, ret$grp)
    ret
}
