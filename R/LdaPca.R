if(FALSE) {

    library(mlrob)
    xx <- get_data("Colon")
    ldap <- LdaPca(xx$x, xx$grp, k=61)
    predict(ldap)
    plot(ldap)

}

LdaPca <- function(x, grp, k=ncol(x), preprocess=c("none", "center", "sphere", "standardize")){

    preprocess <- match.arg(preprocess)

    ## some checks
    clInfo <- class(x)[1]
    if(clInfo == "data.frame") x <- as.matrix(x)

    if(length(grp) != dim(x)[1]){
    stop(paste("grp must be of length", dim(x)[1]))
    }
    if(dim(x)[2] < 1){
    stop("matrix or data.frame expected.")
    }
    n <- nrow(x)
    p <- ncol(x)

    grp <- as.factor(grp)
    glev <- levels(grp)
    ng <- length(glev)


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

##    pca <- PcaClassic(z, k=r)
##    lda <- LdaClassic(getScores(pca), grp)

    pca <- prcomp(z, rank.=r)
    lda <- lda(pca$x, grp)

    pred <- predict(lda)
    grppred <- pred$class

    ## misclassification rate and discriminand scores
    mc <- table(grp, grppred)
    rate <- 1 - sum(diag(mc)) / sum(mc)
    fdiscr <- pred$x               # discriminant scores

    res <- list(pcaobj=pca, ldaobj=lda,
                fdiscr=fdiscr,
                mc=mc, mcrate=rate, grppred=grppred,
                X=x, grp=grp, k=r,
                preprocess=preprocess, scaled_center=scaled_center, scaled_scale=scaled_scale)

    class(res) <- "LdaPca"

  res
}

print.LdaPca <- function(x,...){
  cat("--------------------------------------")
    print(x$ldaobj)

  cat("--------------------------------------\n")
}

plot.LdaPca <- function(obj) {
    proj <- data.frame(obj$pcaobj$x[,1:2])
    proj$grp <- as.factor(obj$grp)
    proj$grppred <- as.factor(obj$grppred)
    firstscores <- NULL
    secondscores <- NULL
    colnames(proj) <- c("firstscores", "secondscores","grp", "grppred")
    gg <- ggplot(proj, aes(firstscores, secondscores, colour = grp, shape = grppred))
    gg <- gg + geom_point()
    gg <- gg + xlab("First PC") + ylab("Second PC")
    gg <- gg + labs(shape="Predicted", color = "Original")
    print(gg)
}

predict.LdaPca <- function(object, newdata){
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
