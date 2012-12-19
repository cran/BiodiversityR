`ensemble.test.nnet` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=5, digits=2, PLOTS=FALSE,    
    Yweights="BIOMOD", factors=NULL, formulae.defaults=TRUE,
    NNET.formula=NULL,
    size=c(2, 4, 6, 8), decay=c(0.1, 0.05, 0.01, 0.001)
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    if (! require(nnet)) {stop("Please install the nnet package")}
    if (formulae.defaults == T) {
        formulae <- ensemble.formulae(x, factors=factors)
    }
    if (is.null(NNET.formula) == T && formulae.defaults == T) {NNET.formula <- formulae$NNET.formula}
    if (is.null(NNET.formula) == T) {stop("Please provide the NNET.formula (hint: use ensemble.formulae function)")}
    environment(NNET.formula) <- .BiodiversityR
    if (is.null(a)==T) {a <- randomPoints(x, n=an, ext=ext)}
    groupp <- kfold(p, k=k)
    groupa <- kfold(a, k=k)
    ns <- length(size)
    nd <- length(decay)
    nt <- ns*nd
    output <- array(NA, dim=c(nt, k+3))
    colnames(output) <- c("size","decay", 1:k,"MEAN")
    output[,"size"] <- rep(size, nd)
    output[,"decay"] <- rep(decay, each=ns) 
    for (i in 1:k){
        cat(paste("\n", "EVALUATION RUN: ", i, "\n", "\n", sep = ""))
        pc <- p[groupp != i,]
        pt <- p[groupp == i,]
        ac <- a[groupa != i,]
        at <- a[groupa == i,]
        assign("pc", pc, envir=.BiodiversityR)
        assign("pt", pt, envir=.BiodiversityR)
        assign("ac", ac, envir=.BiodiversityR)
        assign("at", at, envir=.BiodiversityR)
        for (j in 1:nt) {
            NNET.size <- output[j,"size"]
            NNET.decay <- output[j, "decay"]
            cat(paste("\n", "size: ", NNET.size, ", decay: ", NNET.decay, "\n", sep=""))
            tests <- ensemble.test(x=x, p=pc, a=ac, pt=pt, at=at, 
                PLOTS=PLOTS, evaluations.keep=T,
                MAXENT=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, 
                GAM=0, GAMSTEP=0, MGCV=0, EARTH=0, RPART=0, 
                NNET=1, FDA=0, SVM=0, BIOCLIM=0, DOMAIN=0, MAHAL=0,  
                Yweights=Yweights, factors=factors,
                NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay)
            output[j,2+i] <- tests$NNET.T@auc 
        }
    }
    output[,k+3] <- rowMeans(output[,3:(k+2)], na.rm=T)
    output <- output[order(output[,"MEAN"], decreasing=T),]
    output[,3:(k+3)] <- 100*output[,3:(k+3)]
    output[,3:(k+3)] <- round(output[,3:(k+3)], digits=digits)
    return(output)
}




