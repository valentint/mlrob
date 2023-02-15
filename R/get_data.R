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
##
##  This function uses data sets from the package as well as data from
##  other packages, using the command data(). This will not be liked
##  by "R CMD check" and errors "no visible binding for global variable ..."
##  will be issued. I could not find a solution for this apart from the following:
##  - call data() with a character string, i.e. data("diabetes") instead of data(diabetes)
##  - use the following construct for each data set:
##            data("salmon", package="rrcov", envir=environment)
##            salmon <- get("salmon", envir=environment())
##
get_data <- function(dname=c("Iris", "Wine", "Diabetes", "Crabs",
    "Soil", "Salmon", "Fish", "Olitos", "Fruit", "Ionosphere", "Bupa", "Glass",
    "Pima", "Pima_train", "Pima_test",
    "LetterRecognition", "Thyroid", "Vehicle", "Ecoli", "Heart",
    "MLL", "Gastro", "NIR", "Colon", "Colon_bioconductor", "Golub_bioconductor",
    "Leukemia_big", "Leukemia_small", "Prostate", "SRBCT", "Penicillium",
    "ALL", "ALL_train", "ALL_test")) {

    dname <- match.arg(dname)

    switch(EXPR = dname,
        "Iris" = {
            grp <- datasets::iris$Species
            x <- datasets::iris[, 1:4]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
           cat("\nIris: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Wine" = {
            data("wine27", package="MBCbook", envir=environment())
            wine27 <- get("wine27", envir=environment())
            x <- wine27[, 1:(ncol(wine27)-2)]
            grp <- wine27$Type
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nWine: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Diabetes" = {
            data("diabetes", package="rrcov", envir=environment())
            diabetes <- get("diabetes", envir=environment())
            x <- diabetes[, c("glucose", "insulin", "sspg")]
            grp <- diabetes$group
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nDiabetes (rrcov): n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "MLL" = {
            data("mll")
            x <- mll[, -ncol(mll)]
            grp <- mll$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nMLL: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Gastro" = {
            data("gastro")
            x <- gastro[, -ncol(gastro)]
            grp <- gastro$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nGastro: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "NIR" = {
            data("NIR", package="MBCbook", envir=environment())
           NIR <- get("NIR", envir=environment())
            x <- NIR[, -1]
            grp <- factor(NIR$cls)
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\nNIR: n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },


##
        "Crabs" = {
            data("crabs", package="MASS")
            x <- MASS::crabs[, c("FL", "RW", "CL", "CW", "BD")]
            grp <- MASS::crabs$sex
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Soil" = {
            data("soil", package="rrcov", envir=environment())
            soil <- get("soil", envir=environment())
            x <- soil[soil$D == 0, -c(1,2)]              # only 1983, remove column D (always 0)
            grp <- soil[soil$D == 0, "F"]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Salmon" = {
            data("salmon", package="rrcov", envir=environment())
            salmon <- get("salmon", envir=environment())
            x <- salmon[, -4]
            grp <- salmon$Origin
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Fish" = {
            data("fish", package="rrcov", envir=environment())
            fish <- get("fish", envir=environment())
            # remove observation #14 containing missing value
            myfish <- fish[-14,]
            # The height and width are calculated as percentages
            #   of the third length variable
            myfish[,5] <- myfish[,5] * myfish[,4]/100
            myfish[,6] <- myfish[,6] * myfish[,4]/100
            x <- myfish[, -7]
            grp <- myfish$Species
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Olitos" = {
            data("olitos", package="rrcov", envir=environment())
            olitos <- get("olitos", envir=environment())
            x <- olitos[, -ncol(olitos)]
            grp <- olitos$grp
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Fruit" = {
            data("fruit", package="rrcov", envir=environment())
            fruit <- get("fruit", envir=environment())
            x <- fruit[, -1]
            grp <- fruit$cultivar
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Ionosphere" = {
            data("ionosphere")
            ##  Ionosphere Data Set
            ##  UCI:    https://archive.ics.uci.edu/ml/datasets/ionosphere
            x <- ionosphere[, -ncol(ionosphere)]
            grp <- ionosphere[, ncol(ionosphere)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Bupa" = {
            data("bupa")
            ##  UCI:    https://archive.ics.uci.edu/ml/datasets/liver+disorders
            x <- bupa[, -ncol(bupa)]
            grp <- bupa[, ncol(bupa)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Glass" = {
            ##  Glass Identification Data Set
            ##  UCI:    https://archive.ics.uci.edu/ml/datasets/glass+identification
            ##  mlbench
            data("Glass", package="mlbench", envir=environment())
            Glass <- get("Glass", envir=environment())
            x <- Glass[, -ncol(Glass)]
            grp <- Glass[, ncol(Glass)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Pima" = {
            data("Pima.tr", package="MASS")
            data("Pima.te", package="MASS")
            xx <- rbind(MASS::Pima.tr, MASS::Pima.te)
            x <- xx[, -ncol(xx)]
            grp <- xx[, ncol(xx)]
            rm(xx)
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Pima_train" = {
            data("Pima.tr", package="MASS")
            x <- MASS::Pima.tr[, -ncol(MASS::Pima.tr)]
            grp <- MASS::Pima.tr[, ncol(MASS::Pima.tr)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Pima_test" = {
            data("Pima.te", package="MASS")
            x <- MASS::Pima.te[, -ncol(MASS::Pima.te)]
            grp <- MASS::Pima.te[, ncol(MASS::Pima.te)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "LetterRecognition" = {
            data("LetterRecognition", package="mlbench", envir=environment())
            LetterRecognition <- get("LetterRecognition", envir=environment())
            x <- LetterRecognition[, -1]
            grp <- LetterRecognition[, 1]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Thyroid" = {
            ##  Thyroid Disease Data Set
            ##  UCI:    https://archive.ics.uci.edu/ml/datasets/thyroid+disease
            data("thyroid", package="mclust", envir=environment())
            thyroid <- get("thyroid", envir=environment())
            x <- thyroid[, -1]
            grp <- thyroid[, 1]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Vehicle" = {
            ##  Statlog (Vehicle Silhouettes) Data Set
            ##  UCI:    https://archive.ics.uci.edu/ml/datasets/Statlog+%28Vehicle+Silhouettes%29
            data("vehicle")
            x <- vehicle[, -ncol(vehicle)]
            grp <- vehicle[, ncol(vehicle)]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Ecoli" = {
            ##  Ecoli Data Set
            ##  UCI:    https://archive.ics.uci.edu/ml/datasets/ecoli
            data("ecoli")
            x <- ecoli[, -ncol(ecoli)]
            grp <- ecoli[, ncol(ecoli)]

            ## remove the classes with 5 or less elements
            id <- which(grp %in% names(which(table(grp) <= 5)))
            grp <- grp[-id]
            grp <- factor(grp)
            x <- x[-id,]
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Heart" = {
            ##  Heart Disease Data Set
            ##  UCI:    https://archive.ics.uci.edu/ml/datasets/heart+disease
            data("heartc")
            ## remove the 6 missing values
            heartc <- heartc[!is.na(heartc[,12]) & !is.na(heartc[,13]),]
            x <- heartc[, -ncol(heartc)]
            ## recode class to be 0/1
            grp <- factor(ifelse(heartc[, ncol(heartc)] == 0, 0, 1))
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
##
        "Colon" = {
            data("Colon", package="plsgenomics", envir=environment())
            Colon <- get("Colon", envir=environment())
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

            data("colon")
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

            data("golub_bioconductor")
            x <- as.matrix(golub_bioconductor[, -ncol(golub_bioconductor)])
            grp <- golub_bioconductor$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Leukemia_big" = {
            data("golub")
            x <- as.matrix(leukemia_big[, -ncol(leukemia_big)])
            grp <- leukemia_big$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Leukemia_small" = {
            data("golub")
            x <- as.matrix(leukemia_small[, -ncol(leukemia_small)])
            grp <- leukemia_small$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Prostate" = {
            data("prostate", package="spls", envir=environment())
            prostate <- get("prostate", envir=environment())
            x <- prostate$x
            colnames(x) <- paste0("V", 1:ncol(x))       # caret wants column names of the matrix X
            grp <- as.factor(prostate$y)
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "SRBCT" = {
            data("SRBCT", package="plsgenomics", envir=environment())
            SRBCT <- get("SRBCT", envir=environment())
            x <- SRBCT$X
            colnames(x) <- paste0("V", 1:ncol(x))       # caret wants column names of the matrix X
            grp <- as.factor(SRBCT$Y)
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "Penicillium" = {
            data("penicilliumYES", package="sparseLDA", envir=environment())
            penicilliumYES <- get("penicilliumYES", envir=environment())
            X <- penicilliumYES$X
            x <- X[, -which(apply(X, 2, mad) == 0)]
            grp <- as.factor(c(rep(1,12), rep(2,12), rep(3,12)))
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "ALL_train" = {
            data("Yeoh2002")
            x <- as.matrix(ALL_train[, -ncol(ALL_train)])
            grp <- ALL_train$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "ALL_test" = {
            data("Yeoh2002")
            x <- as.matrix(ALL_test[, -ncol(ALL_test)])
            grp <- ALL_test$Class
            ret <- list(name=dname, x=x, grp=grp, n=nrow(x), p=ncol(x), ng=length(table(grp)))
            cat("\n", dname, ": n=", ret$n, "p=", ret$p, "ng=", ret$ng, "...\n")
        },
        "ALL" = {
            data("Yeoh2002")
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
