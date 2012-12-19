`ensemble.test.gbm` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=5, digits=2, PLOTS=FALSE,    
    Yweights="BIOMOD", factors=NULL,
    GBMSTEP.gbm.x=2:(1+nlayers(x)), 
    complexity=c(3:6), learning=c(0.005, 0.002, 0.001), 
    GBMSTEP.bag.fraction=0.5
)
{
    if (! require(dismo)) {stop("Please install the dismo package")}
    if (! require(gbm)) {stop("Please install the gbm package")}
    if (is.null(a)==T) {a <- randomPoints(x, n=an, ext=ext)}
    groupp <- kfold(p, k=k)
    groupa <- kfold(a, k=k)
    nc <- length(complexity)
    nl <- length(learning)
    nt <- nc*nl
    output <- array(NA, dim=c(nt, 2*k+3))
    colnames(output) <- c("tree.complexity","learning.rate", 1:k,"MEAN",1:k)
    output[,"tree.complexity"] <- rep(complexity, nl)
    output[,"learning.rate"] <- rep(learning, each=nc) 
    for (i in 1:k){
        cat(paste("\n", "EVALUATION RUN: ", i, "\n", "\n", sep = ""))
        pc <- p[groupp != i,]
        pt <- p[groupp == i,]
        ac <- a[groupa != i,]
        at <- a[groupa == i,]
        .BiodiversityR <- new.env()
        assign("pc", pc, envir=.BiodiversityR)
        assign("pt", pt, envir=.BiodiversityR)
        assign("ac", ac, envir=.BiodiversityR)
        assign("at", at, envir=.BiodiversityR)
        for (j in 1:nt) {
            complex <- output[j,"tree.complexity"]
            lr <- output[j, "learning.rate"]
            cat(paste("\n", "complexity: ", complex, ", learning: ", lr, "\n", sep=""))
            tests <- ensemble.test(x=x, p=pc, a=ac, pt=pt, at=at, 
                PLOTS=PLOTS, evaluations.keep=T,
                MAXENT=0, GBM=0, GBMSTEP=1, RF=0, GLM=0, GLMSTEP=0, 
                GAM=0, GAMSTEP=0, MGCV=0, EARTH=0, RPART=0, 
                NNET=0, FDA=0, SVM=0, BIOCLIM=0, DOMAIN=0, MAHAL=0,   
                Yweights=Yweights, factors=factors,   
                GBMSTEP.gbm.x=2:(1+nlayers(x)), 
                GBMSTEP.bag.fraction=GBMSTEP.bag.fraction, 
                GBMSTEP.tree.complexity=complex, GBMSTEP.learning.rate=lr)
            output[j,2+i] <- tests$GBMSTEP.T@auc
            output[j,k+3+i] <- tests$GBMSTEP.trees 
        }
    }
    output[,k+3] <- rowMeans(output[,3:(k+2)], na.rm=T)
    output <- output[order(output[,"MEAN"], decreasing=T),]
    output[,3:(k+3)] <- 100*output[,3:(k+3)]
    output[,3:(k+3)] <- round(output[,3:(k+3)], digits=digits)
    return(output)
}

