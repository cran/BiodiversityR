`check.datasets` <- 
function(x,y) {
    cat("Community dataset:", nrow(x), "X", ncol(x), ", Environmental dataset:", nrow(y), "X", ncol(y), "\n")
    if(nrow(x)!=nrow(y)){
            cat("Warning: community and environmental datasets have different numbers of rows\n")
    }else{
        if(any(rownames(x)!=rownames(y))){
            cat("Warning: rownames for community and environmental datasets are different\n")
        }
    }
}

