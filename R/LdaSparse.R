if(FALSE) {

    library(mlrob)
    ## iris data
    grp <- iris$Species
    x <- iris[, 1:4]

    ldas <- LdaSparse(x, grp)
    predict(ldas)

    xx <- get_data("Gastro")
    ldas <- LdaSparse(xx$x, xx$grp, lambda=0, numVars=14, maxRSS=NA)
    ##  sda(xx$x, xx$grp)
    predict(ldas)

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

    ldas <- LdaSparse(x, grp)
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

LdaSparse <- function(x, grouping, prior=proportions, lambda=1e-6, numVars, maxRSS, maxit=25, trace=FALSE){

    if(is.null(dim(x)))
        stop("x is not a matrix")

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    stopit <- NA
    if((missing(numVars) || is.na(numVars) || is.null(numVars)) &&
       (missing(maxRSS) || is.na(maxRSS) || is.null(maxRSS)))         stopit <- -1     # the default of sda()
    else if(!missing(numVars) && !is.na(numVars) && !is.null(numVars) &&
            !missing(maxRSS) && !is.na(maxRSS) && !is.null(maxRSS))
        stop("Only one of numVars and maxRSS may be provided!")
    else if(!missing(numVars) && !is.na(numVars) && !is.null(numVars)) stopit <- -numVars      # negative for number of variables
    else                                                               stopit <- maxRSS        # positive, to use "penalty"

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

    ##  Normalize the data to mean 0 and unit length columns
    ##      z <- normalize(x)$Xc
    ##  Please note that the scale might be 0 for some variables -
    ##      exclude these variables
    ## Store the result of the normalization in scaled_center and scaled_scale,
    ##  to be used later for predicting of new data.
    ##
    x1 <- scale(x, center=TRUE, scale=FALSE)
    vx <-  sqrt(apply(x1^2, 2, sum))
    id <- which(vx > 0)
    z <- scale(x1[, id], center=FALSE, scale=vx[id])

    scaled_center <- if(!is.null(attr(x1, "scaled:center"))) attr(x1, "scaled:center") else rep(0, p)
    scaled_scale <- vx

    z <- as.data.frame(z)       # to delete the attributes

    ldas <- sda(z, g, lambda=lambda, stop=stopit, maxIte=maxit, trace=trace)
    grppred <- predict(ldas, z)$class

    ## misclassification rate:
    mc <- table(g, grppred)
    ##  mc <- mc[, matchClasses(mc, method = "exact", verbose=FALSE)]
    rate <- 1 - sum(diag(mc)) / sum(mc)

    res <- list(call=xcall, prior=prior, counts=counts,
                lambda=lambda, numVars=ifelse(stopit < 0, -stopit, NA), maxRSS=ifelse(stopit > 0, stopit, NA),
                # meanj=meanj, cv=cv, meanov=meanov,
                mc=mc, mcrate=rate, grppred=grppred,
                fit=ldas,
                method="Sparse LDA",
                X=x, grp=g,
                scaled_center=scaled_center, scaled_scale=scaled_scale)

    class(res) <- "LdaSparse"
    res
}

print.LdaSparse <- function(x,...){
  cat("\n--------------------------------------\n")
  cat("\nResults from Sparse LDA\n")
  print(x$fit)
  cat("\n--------------------------------------\n")
}

plot.LdaSparse <- function(x, ...) {

    #proj <- xc %*%V [,1:2]
    ###proj <- fs[,1:2]

    proj <- data.frame(x$fdiscr)
    proj$grp <- as.factor(x$grp)
    proj$grppred <- as.factor(x$grppred)
    firstscores <- NULL
    secondscores <- NULL
    colnames(proj) <- c("firstscores", "secondscores","grp", "grppred")
    gg <- ggplot(proj, aes(firstscores, secondscores, colour = grp, shape = grppred))
    gg <- gg + geom_point()
    gg <- gg + xlab("First Fisher scores") + ylab("Second Fisher scores")
    print(gg)
}

predict.LdaSparse <- function(object, newdata, ...){
    ct <- FALSE
    if(missing(newdata)) {
        newdata <- object$X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)
    ng <- length(levels(object$grp))    # number of groups

    ## Normalize and predict
    id <- which(object$scaled_scale > 0)
    z <- scale(x[, id, drop=FALSE], center=object$scaled_center[id], scale=object$scaled_scale[id])

    if(ncol(z) != object$fit$origP) stop("Dimensions of training and testing data set different!")

    pred <- predict(object$fit, z)

    ## prediction:
    grp <- pred$class
    grpnam <- as.character(grp)

    ret <- list(grpnam=grpnam, grp=grp, posterior=pred$posterior, x=pred$x)

    if(ct)
        ret$ct <- rrcov::mtxconfusion(object$grp, ret$grp)

    ret
}
