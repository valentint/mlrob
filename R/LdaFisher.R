if(FALSE) {

    ## iris data
    grp <- iris$Species
    x <- iris[, 1:4]

    ldaf <- LdaFisher(x, grp)
    lda <- LdaClassic(x, grp)
    ldap <- LdaPca(x, grp, k=2, preprocess="sphere")
    (tab <- table(grp, predict(ldap, newdata=x)$grpnam)); rrcov:::.AER(tab)


    ## Reclassification results
    (tab <- table(grp, predict(ldaf, newdata=x)$grpnam)); rrcov:::.AER(tab)
    (tab <- table(grp, predict(ldap, newdata=x)$grpnam)); rrcov:::.AER(tab)
    (tab <- table(grp, predict(lda, newdata=x)@classification)); rrcov:::.AER(tab)

    1-cv(ldaf)$aveacc
    1-cv(ldap)$aveacc
    1-cv(lda)$aveacc

    ## diabetes data
    data(diabetes, package="rrcov")
    x <- diabetes[, c("glucose", "insulin", "sspg")]
    grp <- diabetes$group

    ldaf <- LdaFisher(x, grp)
    lda <- LdaClassic(x, grp)
    ldap <- LdaPca(x, grp)

   ## Reclassification results
    (tab <- table(grp, predict(ldaf, newdata=x)$grpnam)); rrcov:::.AER(tab)
    (tab <- table(grp, predict(ldap, newdata=x)$grpnam)); rrcov:::.AER(tab)
    (tab <- table(grp, predict(lda, newdata=x)@classification)); rrcov:::.AER(tab)

    1-cv(ldaf)$aveacc
    1-cv(ldap)$aveacc
    1-cv(lda)$aveacc

    ## wine data
    library(MBCbook)
    data(wine27)

    x <- wine27[, 1:(ncol(wine27)-2)]
    grp <- wine27$Type

    ldaf <- LdaFisher(x, grp)
    lda <- LdaClassic(x, grp)
    ldap <- LdaPca(x, grp, k=2)

   ## Reclassification results
    (tab <- table(grp, predict(ldaf, newdata=x)$grpnam)); rrcov:::.AER(tab)
    (tab <- table(grp, predict(ldap, newdata=x)$grpnam)); rrcov:::.AER(tab)
    (tab <- table(grp, predict(lda, newdata=x)@classification)); rrcov:::.AER(tab)

    1-cv(ldaf)$aveacc
    1-cv(ldap)$aveacc
    1-cv(lda)$aveacc

}

LdaFisher <- function(x, grouping, prior=proportions){

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

    B <- matrix(0, p, p)
    W <- matrix(0, p, p)
    for(j in 1:ng) {
        ###B <- B+prior[j]*((meanj-meanov)%*%t(meanj-meanov))
        B <- B + prior[j] * ((meanj - meanov %*% rep(1, ng)) %*% t(meanj - meanov %*% rep(1, ng)))
        #    W <- W+prior[j]*cov(x[g == lev1[j],])
        W <- W + prior[j] * cv[[j]]
    }

    l <- min(ng-1, p) # use this number of components

    ## Solve the generalized eigenvalue problem using geigen
    V <- geigen(B, W)$vectors
    V <- t(t(V)/(sqrt(diag(t(V) %*% W %*% V))))

    if(FALSE) {
        ## Solve the generalized eigenvalue problem using geigen
        V <- matrix(Re(eigen(solve(W) %*% B)$vec)[, 1:l], ncol=l)
        V <- t(t(V)/(sqrt(diag(t(V) %*% W %*% V))))
    }

    if(FALSE) {
        ## Something from robCompositions, but it seems it does not work -
        ##  see the statement, which is definitely wrong
        ##      l1 <- length(B.svd$d > 1e-6)
        ##  should be:
        ##      l1 <- length(which(B.svd$d > 1e-6))
        ##  but this also does not work

        # besser:
        B.svd <- svd(B)
        l1 <- length(which(B.svd$d > 1e-6))
        B12 <- B.svd$u[, 1:l1] %*% diag(sqrt(B.svd$d[1:l1])) %*% t(B.svd$u[, 1:l1])
        Bm12 <- B.svd$u[, 1:l1] %*% diag(1/sqrt(B.svd$d[1:l1])) %*% t(B.svd$u[, 1:l1])
        K <- eigen(B12 %*% solve(W) %*% t(B12))
        l2 <- min(ng - 1, p)
        l <- min(length(K$val > 1e-6),l2)

        Vs <- Bm12 %*% K$vec[,1:l]
        V <- t(t(Vs) / (sqrt(diag(t(Vs) %*% W %*% Vs))))
    }

    # Fisher scores
    fs <- matrix(NA, nrow=n, ncol=ng)
    dimnames(fs)[[2]] <- lev1
    for(j in 1:ng) {
        xc <- scale(x, meanj[,j], scale=FALSE)
        xproj <- xc %*% V
        fs[,j] <- sqrt(apply(xproj^2, 1, sum) - 2 * log(prior[j]))
    }

    ## prediction:
    grppred <- lev1[apply(fs, 1, which.min)]

    ## misclassification rate:
    mc <- table(g, grppred)
    ##  mc <- mc[, matchClasses(mc, method = "exact", verbose=FALSE)]
    rate <- 1 - sum(diag(mc)) / sum(mc)

    fdiscr <- scale(x, meanov, FALSE) %*% V[, 1:2]               # discriminant scores

    res <- list(call=xcall, prior=prior, counts=counts,
                meanj=meanj, cv=cv, meanov=meanov,
                B=B, W=W, loadings=V,
                scores=fs, fdiscr=fdiscr,
                mc=mc, mcrate=rate, grppred=grppred,
                method="Fisher Linear Discriminant Analysis",
                X=x, grp=g)
    class(res) <- "LdaFisher"
    res
}

print.LdaFisher <- function(x,...){
  cat("--------------------------------------")
  cat("\nResults from Fishers discriminant analysis")
  cat("\n- Variance between the classes: \n")
  print(x$B)
  cat("\n- Variance within the classes: \n")
  print(x$W)
  cat("\n- Loadings matrix: \n")
  print(x$load)
  cat("--------------------------------------\n")
}

plot.LdaFisher <- function(obj) {

    #proj <- xc %*%V [,1:2]
    ###proj <- fs[,1:2]

    proj <- data.frame(obj$fdiscr)
    proj$grp <- as.factor(obj$grp)
    proj$grppred <- as.factor(obj$grppred)
    firstscores <- NULL
    secondscores <- NULL
    colnames(proj) <- c("firstscores", "secondscores","grp", "grppred")
    gg <- ggplot(proj, aes(firstscores, secondscores, colour = grp, shape = grppred))
    gg <- gg + geom_point()
    gg <- gg + xlab("First Fisher scores") + ylab("Second Fisher scores")
    gg <- gg + labs(shape="Predicted", color = "Original")
    print(gg)
}

predict.LdaFisher <- function(object, newdata){
    ct <- FALSE
    if(missing(newdata)) {
        newdata <- object$X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)

    if(length(object$meanj) > 0 & ncol(x) != nrow(object$meanj) | ncol(x) != nrow(object$loadings))
        stop("wrong number of variables")

    ng <- ncol(object$meanj) # number of groups

    # Fisher scores
    fs <- matrix(NA, nrow=nrow(newdata), ncol=ng)
    dimnames(fs)[[2]] <- dimnames(object$mc)[[1]]
    for (j in 1:ng){
        xc <- scale(newdata, object$meanj[,j], scale=FALSE)
        xproj <- xc %*% object$loadings
        fs[,j] <- sqrt(apply(xproj^2, 1, sum) - 2 * log(object$prior[j]))
    }

    ## prediction:
    grp <- apply(fs, 1, which.min)
    grpnam <- colnames(fs)[grp]

    list(grpnam=grpnam, grp=grp)
}
