`diversityresult` <-
function(x,y="",factor,level,index="Shannon",method="all",sortit=F,digits=8){
    jackknife2=function (x, theta, ...) {
        call <- match.call()
        n <- nrow(x)
        u <- rep(0, n)
        for (i in 1:n) {
## changed , in following line from original jackknife of bootstrap library
            u[i] <- theta(x[-i,], ...)
        }
        thetahat <- theta(x, ...)
        jack.bias <- (n - 1) * (mean(u) - thetahat)
        jack.estimate <- thetahat-jack.bias
        jack.se <- sqrt(((n - 1)/n) * sum((u - mean(u))^2))
        return(jack.estimate, jack.se, jack.bias, jack.val = u, call = call)
    }
    diversityresult0=function(x,index="Shannon",method="all"){
        x <- as.matrix(x)
        marg <- 1
        if (index=="Jack.1" || index=="Jack.2" || index=="Chao" || index=="Boot") {method <- "s"}
        if (method=="all") {
            x <- apply(x,2,sum)
            marg <- 2
            }
        INDICES <- c("Shannon","Simpson","inverseSimpson","Logalpha","Berger","richness","abundance","Jevenness",
            "Eevenness","Jack.1","Jack.2","Chao","Boot","ACE")
        index <- match.arg(index, INDICES)
        switch(index, Shannon = {
            result <- diversity(x,index="shannon",MARGIN=marg)    
        }, Simpson = {
            result <- diversity(x,index="simpson",MARGIN=marg)
        }, inverseSimpson = {
            result <- diversity(x,index="invsimpson",MARGIN=marg)
        }, Logalpha = {
            result <- fisher.alpha(x,MARGIN=marg)
        }, Berger = {
            if (marg == 2) {
                result <- max(x)/sum(x)
            }else{
                tots <- as.matrix(apply(x,marg,sum))
                result <- as.matrix(apply(x,marg,max))
                result <- as.matrix(result/tots)[,1]
            }
        }, richness = {
            if (marg == 2) {
                result <- sum(x>0)
            }else{
                result <- as.matrix(apply(x>0,marg,sum))
                result <- result[,1]
            }
        }, abundance = {
            if (marg == 2) {
                result <- sum(x)
            }else{
                result <- as.matrix(apply(x,marg,sum))
                result <- result[,1]
            }
        }, Jevenness = {
            result1 <- diversity(x,index="shannon",MARGIN=marg)
            if (marg == 2) {
                result2 <- sum(x>0)
            }else{
                result2 <- as.matrix(apply(x>0,marg,sum))
                result2 <- result2[,1]
            }
            result <- result1/log(result2)
        }, Eevenness = {
            result1 <- diversity(x,index="shannon",MARGIN=marg)
            if (marg == 2) {
                result2 <- sum(x>0)
            }else{
                result2 <- as.matrix(apply(x>0,marg,sum))
                result2 <- result2[,1]
            }
            result <- exp(result1)/result2
        }, Jack.1 = {
            result2 <- specpool(x)$Jack.1
        }, Jack.2 = {
            result2 <- specpool(x)$Jack.2
        }, Chao = {
            result2 <- specpool(x)$Chao
        }, Boot = {
            result1 <- specpool(x)$Boot
        }, ACE = {
            if (method=="all") {
                result1 <- estimateR(x)["S.ACE"]
            }else{
                result1 <- estimateR(x)["S.ACE",]
            }
        })
    }
    options(digits=digits)   
    if(class(y) == "data.frame") {
        subs <- y[,factor]==level
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q]<-F}
        }
        x <- x[subs,,drop=F]
        freq <- apply(x,2,sum)
        subs <- freq>0
        x <- x[,subs,drop=F]
    }
    x <- as.matrix(x)
    if(dim(x)[1]==0) {
        result <- array(NA,dim=c(1,1))
        colnames(result) <- index
        rownames(result) <- "none"
        return(result)
    }
    if (method=="jackknife") {
        result2 <- jackknife2(x,diversityresult0,index=index)
    }else{
        result <- diversityresult0(x,index=index,method=method)
    }
    if (method=="mean") {
        result2 <- result
        result <- result2[1]
        result[1] <- mean(result2)
    }
    if (method=="sd") {
        result2 <- result
        result <- result2[1]
        result[1] <- sd(result2)
    }
    if (sortit==T && method!="jackknife" && method!="all") {result <- sort(result)}
    if (method!= "jackknife") {
        result2 <- round(result,digits=digits)
        result2 <- data.frame(result2)
        colnames(result2) <- index
    }
    if (method=="all") {rownames(result2) <- "all"}
    if (method=="mean") {rownames(result2) <- "mean"}
    if (method=="sd") {rownames(result2) <- "sd"}
    if (method!="all" && method!="jackknife" && method!="mean" && method!="sd") {rownames(result2) <- names(result)}
    return(result2)
}

