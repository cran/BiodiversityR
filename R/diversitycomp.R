`diversitycomp` <-
function(x,y="",factor1,factor2="",index="Shannon",method="all",sortit=F,...) {
    if (factor2=="") {
        groups <- table(y[,factor1])
        m <- length(groups)
        levels <- names(groups)
        result <- array(NA,dim=c(m,2))
        result[,1] <- groups
        dimnames(result) <- list(factor1=levels,c("n",index))
        names(dimnames(result)) <- c(factor1,"")
        for (i in 1:m) {
            if (method=="all"  || method=="mean"  || method=="sd") {result[i,2] <- diversityresult(x,y,factor1,level=levels[i],method=method,index=index,...)[1,1]}
            if (method=="jackknife") {result[i,2] <- diversityresult(x,y,factor1,level=levels[i],method=method,index=index,...)$jack.estimate}
            if (method!="all" && method!="mean" && method!="sd" && method!="jackknife") {stop(paste("method is not allowed for this function"))}
        }
        if (sortit==T) {
            result2 <- result
            seq <- order(result[,2])
            for (i in 1:m) {
                result[1:m,] <- result2[seq,]
            }
            rownames(result) <- rownames(result2)[seq]
        }
        return(result)
    }else{
        groups <- table(y[,factor1],y[,factor2])
        levels1 <- rownames(groups)
        levels2 <- colnames(groups)
        m1 <- length(levels1)
        m2 <- length(levels2)
        result <- array(NA,dim=c(m1,m2,2))
        result[,,1] <- groups        
        dimnames(result) <- list(factor1=levels1,factor2=levels2,c("n",index))
        names(dimnames(result)) <- c(factor1,factor2,"")
        for (j in 1:m1) {
            subs <- y[,factor1]==levels1[j]
            x1 <- x[subs,,drop=F]
            y1 <- y[subs,,drop=F]
            for (i in 1:m2) {
                if (method=="all" || method=="mean"  || method=="sd") {result[j,i,2] <- diversityresult(x1,y1,factor2,level=levels2[i],method=method,index=index,...)[1,1]}
                if (method=="jackknife") {result[j,i,2] <- diversityresult(x1,y1,factor2,level=levels2[i],method=method,index=index,...)$jack.estimate}
            }
        }
        return(result)
    }
}

