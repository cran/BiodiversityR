`check.datasets` <- 
function(x,y) {
    nfact <- 0
    factors <- NULL
    all.good <- TRUE
    for (i in 1:ncol(x)) {
        if(is.factor(x[,i])) {
            nfact <- nfact+1
            factors <- c(factors, colnames(x)[i])
        }
    }
    if (nfact>0) {
        all.good <- FALSE
        cat("Warning:", nfact, "variables of the community dataset ( out of a total of", ncol(x), ") are factors\n")
        print(factors)
    }
    if(nrow(x)!=nrow(y)){
        all.good <- FALSE
        cat("Warning: community and environmental datasets have different numbers of rows\n")
    }else{
        if(any(rownames(x)!=rownames(y))){
            all.good <- FALSE
            cat("Warning: rownames for community and environmental datasets are different\n")
        }
    }
    if (all.good == TRUE) {cat("OK\n")}
}

