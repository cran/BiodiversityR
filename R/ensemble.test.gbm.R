`ensemble.test.gbm` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, excludep=FALSE, ext=NULL, k=5, 
    TrainData=NULL,
    layer.drops=NULL, VIF=FALSE,
    PLOTS=FALSE,    
    Yweights="BIOMOD", factors=NULL,
    GBMSTEP.gbm.x=2:(ncol(TrainData.orig)), 
    complexity=c(3:6), learning=c(0.005, 0.002, 0.001), 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    if (! require(gbm)) {stop("Please install the gbm package")}
    k <- as.integer(k)
    if (k < 2) {
        cat(paste("\n", "NOTE: parameter k was set to be smaller than 2", sep = ""))
        cat(paste("\n", "default value of 5 therefore set for parameter k", "\n", sep = ""))
        k <- 5
    }
# check data
    if (is.null(TrainData) == T) {
        if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
        if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
        if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
    }
# check TrainData
    if (is.null(TrainData) == F) {
        TrainData <- data.frame(TrainData)
        if (colnames(TrainData)[1] !="pb") {stop("first column for TrainData should be 'pb' containing presence (1) and absence (0) data")}
        if ((is.null(x) == F) && (nlayers(x) != (ncol(TrainData)-1))) {
            cat(paste("\n", "WARNING: different number of explanatory variables in rasterStack and data.frame", sep = ""))
        }
    }
# modify rasterStack x only if TrainData is not used
    if (is.null(TrainData) == T) {
        if (is.null(layer.drops) == F) {
            vars <- names(x)
            layer.drops <- as.character(layer.drops)
            factors <- as.character(factors)
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                if (any(vars==layer.drops[i])==FALSE) {
                    cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among grid layers", "\n", sep = ""))
                }else{
                    cat(paste("\n", "NOTE: variable ", layer.drops[i], " will not be included as explanatory variable", "\n", sep = ""))
                    x <- dropLayer(x, which(names(x) == layer.drops[i]))
                    vars <- names(x)
                    if (is.null(factors) == F) {
                        factors <- factors[factors != layer.drops[i]]
                        if(length(factors) == 0) {factors <- NULL}
                    }
                }
            }
        }
        if (is.null(factors) == F) {
            vars <- names(x)
            factors <- as.character(factors)
            nf <- length(factors)
            for (i in 1:nf) {
                if (any(vars==factors[i])==FALSE) {
                    cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among grid layers", "\n", sep = ""))
                }
            }
        }
        # set minimum and maximum values
            for (i in 1:nlayers(x)) {
                x[[i]] <- setMinMax(x[[i]])
            }
        # declare factor layers
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                j <- which(names(x) == factors[i])
                x[[j]] <- raster::as.factor(x[[j]])
            }
        }
# modify TrainData
    }else{
        if (is.null(layer.drops) == F) {
            vars <- colnames(TrainData)
            layer.drops <- as.character(layer.drops)
            factors <- as.character(factors)
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                if (any(vars==layer.drops[i])==FALSE) {
                    cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among columns of TrainData", "\n", sep = ""))
                }else{
                    cat(paste("\n", "NOTE: variable ", layer.drops[i], " will not be included as explanatory variable", "\n", sep = ""))
                    TrainData <- TrainData[, which(colnames(TrainData) != layer.drops[i])]
                    vars <- colnames(TrainData)
                    if (is.null(factors) == F) {
                        factors <- factors[factors != layer.drops[i]]
                        if(length(factors) == 0) {factors <- NULL}
                    }
                }
            }
        }
        if (is.null(factors) == F) {
            vars <- colnames(TrainData)
            factors <- as.character(factors)
            nf <- length(factors)
            for (i in 1:nf) {
                if (any(vars==factors[i])==FALSE) {
                     cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among columns of TrainData", "\n", sep = ""))
                }
            }
        }
    }

# create TrainData and TestData
    if (is.null(TrainData) == F) {
        if(any(is.na(TrainData))) {
            cat(paste("\n", "WARNING: sample units with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData)
        TrainData <- TrainData[TrainValid,]
    }else{
        if (is.null(a)==T) {
            if (excludep == T) {
                a <- randomPoints(x[[1]], n=an, p=p, ext=ext, excludep=T)
            }else{
                a <- randomPoints(x[[1]], n=an, p=NULL, ext=ext, excludep=F)
            }        
        }
        TrainData <- prepareData(x, p, b=a, factors=factors, xy=FALSE)
        if(any(is.na(TrainData[TrainData[,"pb"]==1,]))) {
            cat(paste("\n", "WARNING: presence locations with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==1,])
        p <- p[TrainValid,]
        if(any(is.na(TrainData[TrainData[,"pb"]==0,]))) {
            cat(paste("\n", "WARNING: background locations with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==0,])
        a <- a[TrainValid,]
        TrainData <- prepareData(x, p, b=a, factors=factors, xy=FALSE)
    }
    TrainData.orig <- TrainData
    assign("TrainData.orig", TrainData.orig, envir=.BiodiversityR)
#
    nc <- length(complexity)
    nl <- length(learning)
    nt <- nc*nl
    output <- array(NA, dim=c(nt, 2*k+3))
    colnames(output) <- c("tree.complexity","learning.rate", 1:k,"MEAN",1:k)
    output[,"tree.complexity"] <- rep(complexity, nl)
    output[,"learning.rate"] <- rep(learning, each=nc) 
#
    groupp <- kfold(TrainData, k=k, by=TrainData[,"pb"])
    for (i in 1:k){
        cat(paste("\n", "EVALUATION RUN: ", i, "\n\n", sep = ""))
        TrainData.c <- TrainData[groupp != i,]
        TestData.c <- TrainData[groupp == i,]
        assign("TrainData.c", TrainData.c, envir=.BiodiversityR)
        assign("TestData.c", TestData.c, envir=.BiodiversityR)
        for (j in 1:nt) {
            complex <- output[j,"tree.complexity"]
            lr <- output[j, "learning.rate"]
            cat(paste("\n", "complexity: ", complex, ", learning: ", lr, "\n", sep=""))
            tests <- ensemble.test(x=x,
                TrainData=TrainData.c, TestData=TestData.c, 
                VIF=VIF,
                PLOTS=PLOTS, evaluations.keep=T,
                MAXENT=0, GBM=0, GBMSTEP=1, RF=0, GLM=0, GLMSTEP=0, 
                GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, 
                NNET=0, FDA=0, SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, GEODIST=0,   
                Yweights=Yweights, factors=factors,   
                GBMSTEP.gbm.x=2:(1+nlayers(x)), 
                GBMSTEP.bag.fraction=GBMSTEP.bag.fraction, 
                GBMSTEP.tree.complexity=complex, GBMSTEP.learning.rate=lr, 
                GBMSTEP.step.size=GBMSTEP.step.size)
            output[j,2+i] <- tests$GBMSTEP.T@auc
            output[j,k+3+i] <- tests$GBMSTEP.trees 
        }
    }
    output[,k+3] <- rowMeans(output[,3:(k+2)], na.rm=T)
    output <- output[order(output[,"MEAN"], decreasing=T),]
    output[,3:(k+3)] <- 100*output[,3:(k+3)]
    print(output)
    cat(paste("\n\n"))
    return(list(table=output, call=match.call() ))
}


