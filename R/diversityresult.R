`diversityresult` <- function(
    x, y=NULL, factor=NULL, level=NULL, 
    index=c("Shannon", "Simpson", "inverseSimpson", "Logalpha", "Berger", 
        "richness", "abundance", "Jevenness", "Eevenness", 
        "jack1", "jack2", "chao", "boot"),
    method=c("pooled", "each site", "mean", "sd", "jackknife"),
    sortit=FALSE, digits=8)
{

    INDEX <- c("Shannon", "Simpson", "inverseSimpson", "Logalpha", "Berger", 
        "richness", "abundance", "Jevenness", "Eevenness", 
        "jack1", "jack2", "chao", "boot")
    if ((index %in% INDEX) == F) {stop(paste("choose an accepted index, not index: ", index, sep=""))}

    METHOD <- c("pooled", "each site", "mean", "sd", "jackknife")
    if ((method %in% METHOD) == F) {stop(paste("choose an accepted method, not method: ", method, sep=""))}

    if (is.null(y) == F) {
        if((factor %in% names(y)) == F) {stop("specified factor '", factor, "' is not a variable of the environmental data frame")}
        if(is.factor(y[, factor]) == F) {stop("specified factor '", factor, "' is not a factor")} 
        levels1 <- as.character(levels(as.factor(as.character(y[, factor]))))
        if((level %in% levels1) == F) {stop("specified level '", level, "' is not an available factor level")}
    }

    if (index %in% c("jack1", "jack2", "chao", "boot") && method != "pooled") {
        cat(paste("\n", "Note that default method for index '", index, "' is method 'pooled'", "\n\n", sep=""))
        method <- "pooled"
    }

    diversityresult0=function(
        x, index, method)
    {
        marg <- 1
        if (method=="pooled" && index!="jack1" && index!="jack2" && index!="chao" && index!="boot") {
            x <- apply(x, 2, sum)
            marg <- 2
        }
        if (index == "Shannon") {result <- diversity(x, index="shannon", MARGIN=marg)}
        if (index == "Simpson") {result <- diversity(x, index="simpson", MARGIN=marg)}
        if (index == "inverseSimpson") {result <- diversity(x, index="invsimpson", MARGIN=marg)}
        if (index == "Logalpha") {result <- fisher.alpha(x, MARGIN=1)}
        if (index == "Berger") {
            if (marg == 2) {
                result <- max(x)/sum(x)
            }else{
                tots <- as.matrix(apply(x, marg, sum))
                result <- as.matrix(apply(x, marg, max))
                result <- as.matrix(result/tots)[,1]
            }
        }
        if (index == "richness") {
            if (marg == 2) {
                result <- sum(x>0)
            }else{
                result <- as.matrix(apply(x>0, marg, sum))
                result <- result[,1]
            }
        }
         if (index == "abundance") {
            if (marg == 2) {
                result <- sum(x)
            }else{
                result <- as.matrix(apply(x, marg, sum))
                result <- result[,1]
            }
        }
        if (index == "Jevenness") {
            result1 <- diversity(x, index="shannon", MARGIN=marg)
            if (marg == 2) {
                result2 <- sum(x>0)
            }else{
                result2 <- as.matrix(apply(x>0, marg, sum))
                result2 <- result2[,1]
            }
            result <- result1/log(result2)
        }

        if (index == "Eevenness") {
            result1 <- diversity(x, index="shannon", MARGIN=marg)
            if (marg == 2) {
                result2 <- sum(x>0)
            }else{
                result2 <- as.matrix(apply(x>0, marg, sum))
                result2 <- result2[,1]
            }
            result <- exp(result1)/result2
        }
        if (index == "jack1") {result <- specpool(x)$jack1}
        if (index == "jack2") {result <- specpool(x)$jack2}
        if (index == "chao") {result <- specpool(x)$chao}
        if (index == "boot") {result <- specpool(x)$boot}

        return(result)
    }

    options(digits=digits)   

    if(is.null(y) == F) {
        subs <- y[, factor] == as.character(level)
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q]<-F}
        }
        x <- x[subs, , drop=F]
        freq <- apply(x, 2, sum)
        subs <- freq > 0
        x <- x[, subs, drop=F]
    }

    x <- as.matrix(x)
    if(dim(x)[1]==0) {
        result <- array(NA,dim=c(1,1))
        colnames(result) <- index
        rownames(result) <- "none"
        return(result)
    }
    if (method == "jackknife") {
#        if (! require(bootstrap)) {stop("Please install the bootstrap package")}
        thetadiv <- function(x, xdata, index) {
            xdata2 <- xdata[x, 1:ncol(xdata)] 
            diversityresult0(xdata2, index=index, method="pooled")
        }
        if (nrow(x) > 1) {
            result2 <- bootstrap::jackknife(1:nrow(x), thetadiv, x, index=index)
            result2$jack.estimate <- mean(as.numeric(result2$jack.values), na.rm=T)
        }else{
            result2 <- list(jack.values=NA, jack.estimate=NA)
        }
    }
    if (method != "jackknife") {
        result <- diversityresult0(x, index=index, method=method)
    }
    if (method == "mean") {
        result2 <- result
        result <- result2[1]
        result[1] <- mean(result2)
    }
    if (method == "sd") {
        result2 <- result
        result <- result2[1]
        result[1] <- sd(result2)
    }
    if (sortit == T && method != "jackknife" && method != "pooled") {result <- sort(result)}
    if (method!="jackknife") {
        result2 <- round(result, digits=digits)
        result2 <- data.frame(result2)
        colnames(result2) <- index
    }
    if (method=="pooled") {rownames(result2) <- "pooled"}
    if (method=="mean") {rownames(result2) <- "mean"}
    if (method=="sd") {rownames(result2) <- "sd"}
    if (method!="pooled" && method!="jackknife" && method!="mean" && method!="sd") {rownames(result2) <- names(result)}

    return(result2)
}

