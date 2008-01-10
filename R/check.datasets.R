`check.datasets` <- 
function(x,y) {
    nfact <- 0
	for (i in 1:ncol(x)){
        if(is.factor(x[,i])) {nfact <- nfact+1}
	}
    if (nfact>0){
        cat("Warning:", nfact, "variables of the community dataset ( out of a total of", ncol(x), ") are factors\n")
    }
    if(nrow(x)!=nrow(y)){
        cat("Warning: community and environmental datasets have different numbers of rows\n")
    }else{
        if(any(rownames(x)!=rownames(y))){
            cat("Warning: rownames for community and environmental datasets are different\n")
        }
    }
}

