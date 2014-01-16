`ensemble.strategy` <- function(
    TrainData=NULL, 
    strategy.k=4, strategy.runs=10, verbose=FALSE,
    ENSEMBLE.best=c(1:19), ENSEMBLE.min=c(0.7),
    ENSEMBLE.decay=c(1.0, 1.25, 1.5, 1.75, 2.0), ENSEMBLE.interval.width=c(0.05) 
)
{
    if (! require(dismo)) {stop("Please install the dismo package")}
#   input AUC
    modelnames <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
    weights <- numeric(length=length(modelnames))
    final.weights <- weights
    names(weights) <- modelnames
    bests <- length(ENSEMBLE.best)
    decays <- length(ENSEMBLE.decay)
    intervals <- length(ENSEMBLE.interval.width)
    mins <- length(ENSEMBLE.min) 
#
# output for each cross-validation run
    output <- data.frame(array(dim=c(bests*decays*intervals*mins, 8), NA))
    if (nrow(output) == 1) {cat(paste("\n", "NOTE: no alternatives available for choosing best strategy", "\n", sep=""))}
    colnames(output) <- c("ENSEMBLE.best", "ENSEMBLE.decay", "ENSEMBLE.interval.width", "ENSEMBLE.min", "model.C", "AUC.C", "model.T", "AUC.T")
    all.combinations <- expand.grid(ENSEMBLE.best, ENSEMBLE.decay, ENSEMBLE.interval.width, ENSEMBLE.min)
    output[,c(1:4)] <- all.combinations
#
# repeat for the internal cross-validation runs
    for (s in 1:strategy.runs) {
    groupd <- kfold(TrainData, k=strategy.k, by=TrainData[,"pb"])
    weights.s <- weights

# k-fold cross-validation
    for (c in 1:strategy.k) { 

    TrainData1 <- TrainData[groupd != c,]
    TestData1 <- TrainData[groupd == c,]

    weights.cal <- c(0, weights)
    weights.eval <- c(0, weights)
    names(weights.cal) <- c("ENSEMBLE", modelnames)
    names(weights.eval) <- c("ENSEMBLE", modelnames)
    for (i in 1:length(weights)) {
        TrainPres <- TrainData1[TrainData1[,"pb"]==1, modelnames[i]]
        TrainAbs <- TrainData1[TrainData1[,"pb"]==0, modelnames[i]]
        if (sum(TrainPres, TrainAbs, na.rm=T) != 0) {
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            weights.cal[i+1] <- eval1@auc
        }
        TestPres <- TestData1[TestData1[,"pb"]==1, modelnames[i]]
        TestAbs <- TestData1[TestData1[,"pb"]==0, modelnames[i]]
        if (sum(TestPres, TestAbs, na.rm=T) != 0) {
            eval2 <- evaluate(p=TestPres, a=TestAbs)
            weights.eval[i+1] <- eval2@auc
        }
    }
    input.weights.c <- weights.cal[names(weights.cal) != "ENSEMBLE"]
#    if (verbose == T) {
#        cat(paste("\n", "input calibration weights", "\n", sep=""))
#        print(input.weights.c)
#    }
    input.weights.e <- weights.eval[names(weights.eval) != "ENSEMBLE"]
#    if (verbose == T) {
#        cat(paste("\n", "input evaluation weights", "\n", sep=""))
#        print(input.weights.e)
#    }
    auc.target <- -1.0
    for (r in 1:nrow(output)) {
#        if (verbose == T) {
#            cat(paste("\n", "run ", r, ": best=", output[r, "ENSEMBLE.best"], ": decay=", output[r, "ENSEMBLE.decay"],  ", interval.width=", output[r, "ENSEMBLE.interval.width"], ", min=", output[r, "ENSEMBLE.min"] ,"\n", sep=""))
#        }
#
# strategy based on evaluations
        ws <- ensemble.weights(input.weights.e, decay=output[r, "ENSEMBLE.decay"], best=output[r, "ENSEMBLE.best"], 
            interval.width=output[r, "ENSEMBLE.interval.width"], min.weight=output[r, "ENSEMBLE.min"])
#        if (verbose == T) {print(ws)}
        TrainData1[,"ENSEMBLE"] <- ws["MAXENT"]*TrainData1[,"MAXENT"] + ws["GBM"]*TrainData1[,"GBM"] +
            ws["GBMSTEP"]*TrainData1[,"GBMSTEP"] + ws["RF"]*TrainData1[,"RF"] + ws["GLM"]*TrainData1[,"GLM"] +
            ws["GLMSTEP"]*TrainData1[,"GLMSTEP"] + ws["GAM"]*TrainData1[,"GAM"] + ws["GAMSTEP"]*TrainData1[,"GAMSTEP"] +
            ws["MGCV"]*TrainData1[,"MGCV"] + ws["MGCVFIX"]*TrainData1[,"MGCVFIX"] + ws["EARTH"]*TrainData1[,"EARTH"] +
            ws["RPART"]*TrainData1[,"RPART"] + ws["NNET"]*TrainData1[,"NNET"] + ws["FDA"]*TrainData1[,"FDA"] +
            ws["SVM"]*TrainData1[,"SVM"] + ws["SVME"]*TrainData1[,"SVME"] + ws["BIOCLIM"]*TrainData1[,"BIOCLIM"] +
            ws["DOMAIN"]*TrainData1[,"DOMAIN"] + ws["MAHAL"]*TrainData1[,"MAHAL"]
#        TrainData[,"ENSEMBLE"] <- trunc(TrainData[,"ENSEMBLE"])
        TrainPres <- TrainData1[TrainData1[,"pb"]==1,"ENSEMBLE"]
        TrainAbs <- TrainData1[TrainData1[,"pb"]==0,"ENSEMBLE"]
        eval1 <- NULL
        if (sum(TrainPres, TrainAbs, na.rm=T) != 0) {
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
#            if (verbose == T) {
#                cat(paste("\n", "evaluation with train data", "\n", sep=""))
#                print(eval1)
#            }
            weights.cal["ENSEMBLE"] <- eval1@auc
            weights.cal <- weights.cal[order(weights.cal, decreasing=T)]
            output[r, "model.C"] <- names(weights.cal)[1]
            output[r, "AUC.C"] <- weights.cal[1]
        }
        TestData1[,"ENSEMBLE"] <- ws["MAXENT"]*TestData1[,"MAXENT"] + ws["GBM"]*TestData1[,"GBM"] +
            ws["GBMSTEP"]*TestData1[,"GBMSTEP"] + ws["RF"]*TestData1[,"RF"] + ws["GLM"]*TestData1[,"GLM"] +
            ws["GLMSTEP"]*TestData1[,"GLMSTEP"] + ws["GAM"]*TestData1[,"GAM"] + ws["GAMSTEP"]*TestData1[,"GAMSTEP"] +
            ws["MGCV"]*TestData1[,"MGCV"] + ws["MGCVFIX"]*TestData1[,"MGCVFIX"] + ws["EARTH"]*TestData1[,"EARTH"] +
            ws["RPART"]*TestData1[,"RPART"] + ws["NNET"]*TestData1[,"NNET"] + ws["FDA"]*TestData1[,"FDA"] +
            ws["SVM"]*TestData1[,"SVM"] + ws["SVME"]*TestData1[,"SVME"] + ws["BIOCLIM"]*TestData1[,"BIOCLIM"] +
            ws["DOMAIN"]*TestData1[,"DOMAIN"] + ws["MAHAL"]*TestData1[,"MAHAL"]
#        TestData[,"ENSEMBLE"] <- trunc(TestData[,"ENSEMBLE"])
        TestPres <- TestData1[TestData1[,"pb"]==1,"ENSEMBLE"]
        TestAbs <- TestData1[TestData1[,"pb"]==0,"ENSEMBLE"]
        eval2 <- NULL
        if (sum(TestPres, TestAbs, na.rm=T) != 0) {
            eval2 <- evaluate(p=TestPres, a=TestAbs)
#            if (verbose == T) {
#                cat(paste("\n", "evaluation with test data", "\n", sep=""))
#                print(eval2)
#            }
            weights.eval["ENSEMBLE"] <- eval2@auc
            weights.eval <- weights.eval[order(weights.eval, decreasing=T)]
            output[r, "model.T"] <- names(weights.eval)[1]
            output[r, "AUC.T"] <- weights.eval[1]
        }
        if (weights.eval[1] > auc.target) {
            auc.target <- weights.eval[1]
            weights.out <- ws
        }
    }
    output <- output[order(output[,"AUC.T"], decreasing=T), ]
    if (verbose==T && is.na(output[1, "AUC.T"])==F) {
        cat(paste("\n", "RUN: ", s, "; CROSSVALIDATION: ", c, "; Best strategy", "\n", sep = ""))
        print(as.vector(output[1,]))
    }
    if (verbose == T) {
        cat(paste("\n", "Weights used for best strategy", "\n", sep = ""))
        print(weights.out)
    }

    weights.s <- weights.s + weights.out/strategy.k

# end of k-fold
}

    if (verbose == T) {
        cat(paste("\n", "RUN: ", s, ": Average of weights used for best strategy", "\n", sep = ""))
        print(weights.s)
    }

    final.weights <- final.weights + weights.s/strategy.runs

# end of runs
}

    cat(paste("\n", "Average of weights used for best strategy (internal cross-validations: ", strategy.k, ", runs: ", strategy.runs, ") \n", sep = ""))
    print(final.weights)
    return(list(weights=final.weights))
}

