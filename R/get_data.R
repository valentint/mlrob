if(FALSE) {

    dd <- get_data("Ionosphere1")
    ldar <- LdaRotation(dd$x, dd$grp, preprocess="sphere")
    ldar$mc

    (tab <- table(dd$grp, predict(ldar, newdata=dd$x)$grp)); rrcov:::.AER(tab)
    1-cv(ldar)$aveacc
    loo(ldar)$eaer
    holdout(ldar)

    ldaf <- LdaFisher(dd$x, dd$grp)

    (tab <- table(dd$grp, predict(ldaf, newdata=dd$x)$grp)); rrcov:::.AER(tab)
    1-cv(ldaf)$aveacc
    loo(ldaf)$eaer
    holdout(ldaf)

}

get_data <- function(dname=c("Iris", "Wine", "Diabetes", "Crabs",
    "Soil", "Salmon", "Fish", "Olitos", "Fruit", "Ionosphere", "Bupa", "Glass",
    "Pima", "Pima_train", "Pima_test",
    "LetterRecognition", "Thyroid",
    "MLL", "Gastro", "NIR", "Colon", "Colon_bioconductor", "Golub_bioconductor",
    "Leukemia_big", "Leukemia_small", "Prostate", "SRBCT", "Penicillium",
    "ALL", "ALL_train", "ALL_test")) {

    dname <- match.arg(dname)

    switch(EXPR = dname,
        "Iris" = {
            grp <- iris$Species
            x <- iris[, 1:4]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
           cat("\nIris: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Wine" = {
            data(wine27, package="MBCbook")

            x <- wine27[, 1:(ncol(wine27)-2)]
            grp <- wine27$Type
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nWine: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Diabetes" = {
            data(diabetes, package="rrcov")
            x <- diabetes[, c("glucose", "insulin", "sspg")]
            grp <- diabetes$group
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nDiabetes (rrcov): n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "MLL" = {
            data(mll)
            x <- mll[, -ncol(mll)]
            grp <- mll$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nMLL: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Gastro" = {
            data(gastro)
            x <- gastro[, -ncol(gastro)]
            grp <- gastro$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nGastro: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "NIR" = {
            data(NIR, package="MBCbook")
            x <- NIR[, -1]
            grp <- factor(NIR$cls)
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nNIR: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },


##
        "Crabs" = {
            data(crabs, package="MASS")
            x <- crabs[, c("FL", "RW", "CL", "CW", "BD")]
            grp <- crabs$sex
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Soil" = {
            data(soil, package="rrcov")
            x <- soil[soil$D == 0, -c(1,2)]              # only 1983, remove column D (always 0)
            grp <- soil[soil$D == 0, "F"]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Salmon" = {
            data(salmon, package="rrcov")
            x <- salmon[, -4]
            grp <- salmon$Origin
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Fish" = {
            data(fish, package="rrcov")
            # remove observation #14 containing missing value
            fish <- fish[-14,]
            # The height and width are calculated as percentages
            #   of the third length variable
            fish[,5] <- fish[,5]*fish[,4]/100
            fish[,6] <- fish[,6]*fish[,4]/100
            x <- fish[, -7]
            grp <- fish$Species
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Olitos" = {
            data(olitos, package="rrcov")
            x <- olitos[, -ncol(olitos)]
            grp <- olitos$grp
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Fruit" = {
            data(fruit, package="rrcov")
            x <- fruit[, -1]
            grp <- fruit$cultivar
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Ionosphere" = {
            data(ionosphere)
            x <- ionosphere[, -ncol(ionosphere)]
            grp <- ionosphere[, ncol(ionosphere)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Bupa" = {
            data(bupa)
            x <- bupa[, -ncol(bupa)]
            grp <- bupa[, ncol(bupa)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Glass" = {
            data(Glass, package="mlbench")
            x <- Glass[, -ncol(Glass)]
            grp <- Glass[, ncol(Glass)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Pima" = {
            data(Pima.tr, package="MASS")
            data(Pima.te, package="MASS")
            xx <- rbind(Pima.tr, Pima.te)
            x <- xx[, -ncol(xx)]
            grp <- xx[, ncol(xx)]
            rm(xx)
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Pima_train" = {
            data(Pima.tr, package="MASS")
            x <- Pima.tr[, -ncol(Pima.tr)]
            grp <- Pima.tr[, ncol(Pima.tr)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Pima_test" = {
            data(Pima.te, package="MASS")
            x <- Pima.te[, -ncol(Pima.te)]
            grp <- Pima.te[, ncol(Pima.te)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "LetterRecognition" = {
            data(LetterRecognition, package="mlbench")
            x <- LetterRecognition[, -1]
            grp <- LetterRecognition[, 1]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Thyroid" = {
            data(thyroid, package="mclust")
            x <- thyroid[, -1]
            grp <- thyroid[, 1]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
##
        "Colon" = {
            data(Colon, package="plsgenomics")
            x <- Colon$X
            grp <- as.factor(Colon$Y)
##            if(prep) {
##                x <- log(x)
##                xmed <- apply(x, 2, median)
##                x <- sweep(x, 2, xmed)
##            }
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Colon_bioconductor" = {
            ##library(Biobase)
            ##library(colonCA)
            ##data(colonCA)
            ##x <- t(exprs(colonCA))
            ##grp <- colonCA$class

            data(colon)
            x <- colon[, -ncol(colon)]
            grp <- colon$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Golub_bioconductor" = {
            ##library(Biobase)
            ##library(golubEsets)
            ##data(Golub_Merge)
            ##x <- t(exprs(Golub_Merge))
            ##grp <- Golub_Merge$ALL.AML

            data(golub_bioconductor)
            x <- as.matrix(golub_bioconductor[, -ncol(golub_bioconductor)])
            grp <- golub_bioconductor$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Leukemia_big" = {
            data(golub)
            x <- as.matrix(leukemia_big[, -ncol(leukemia_big)])
            grp <- leukemia_big$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Leukemia_small" = {
            data(golub)
            x <- as.matrix(leukemia_small[, -ncol(leukemia_small)])
            grp <- leukemia_small$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Prostate" = {
            data(prostate, package="spls")
            x <- prostate$x
            colnames(x) <- paste0("V", 1:ncol(x))       # caret wants column names of the matrix X
            grp <- as.factor(prostate$y)
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "SRBCT" = {
            data(SRBCT, package="plsgenomics")
            x <- SRBCT$X
            colnames(x) <- paste0("V", 1:ncol(x))       # caret wants column names of the matrix X
            grp <- as.factor(SRBCT$Y)
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Penicillium" = {
            data(penicilliumYES, package="sparseLDA")
            X <- penicilliumYES$X
            x <- X[, -which(apply(X, 2, mad) == 0)]
            grp <- as.factor(c(rep(1,12), rep(2,12), rep(3,12)))
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "ALL_train" = {
            data(Yeoh2002)
            x <- as.matrix(ALL_train[, -ncol(ALL_train)])
            grp <- ALL_train$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "ALL_test" = {
            data(Yeoh2002)
            x <- as.matrix(ALL_test[, -ncol(ALL_test)])
            grp <- ALL_test$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "ALL" = {
            data(Yeoh2002)
            xx <- rbind(ALL_train, ALL_test)
            x <- as.matrix(xx[, -ncol(xx)])
            grp <- xx$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        {
            cat("\nData set ", dname, "not found! \n")
            ret <- NULL
        }

        )
    return(ret)
}
