`diversitycomp` <- function(
    x, y=NULL, factor1=NULL, factor2=NULL, 
    index=c("Shannon", "Simpson", "inverseSimpson", "Logalpha", "Berger",
        "simpson.unb", "simpson.unb.inverse", 
        "richness", "abundance", "Jevenness", "Eevenness", 
        "jack1", "jack2", "chao", "boot"),
    method=c("pooled", "mean", "sd", "max", "jackknife"), 
    sortit=FALSE, digits=8) 
{

    INDEX <- c("Shannon", "Simpson", "inverseSimpson", "Logalpha", "Berger", 
        "simpson.unb", "simpson.unb.inverse",
        "richness", "abundance", "Jevenness", "Eevenness", 
        "jack1", "jack2", "chao", "boot")
    if ((index %in% INDEX) == F) {stop(paste("choose an accepted index, not index: ", index, sep=""))}

    METHOD <- c("pooled", "mean", "sd", "max", "jackknife")
    if ((method %in% METHOD) == F) {stop(paste("choose an accepted method, not method: ", method, sep=""))}

    if (is.null(y) == F) {
        if((factor1 %in% names(y)) == F) {stop("specified factor1 '", factor1, "' is not a variable of the environmental data frame")}
        if(is.factor(y[, factor1]) == F) {stop("specified factor1 '", factor1, "' is not a factor")}
        y[, factor1] <- as.factor(as.character(y[, factor1]))
        if (is.null(factor2) == F) {
            if((factor2 %in% names(y)) == F) {stop("specified factor2 '", factor2, "' is not a variable of the environmental data frame")}
            if(is.factor(y[, factor2]) == F) {stop("specified factor2 '", factor2, "' is not a factor")} 
            y[, factor2] <- as.factor(as.character(y[, factor2]))
           }
    }

    if (is.null(factor2) == T) {
        groups <- table(y[, factor1])
        m <- length(groups)
        levels <- as.character(names(groups))
        result <- array(NA, dim=c(m, 2))
        result[, 1] <- groups
        dimnames(result) <- list(factor1=levels, c("n", index))
        names(dimnames(result)) <- c(factor1, "")
        for (i in 1:m) {
            if (method %in% c("pooled", "mean", "max", "sd")) {result[i, 2] <- as.numeric(diversityresult(x, y, factor=factor1, level=levels[i], method=method, index=index, digits=digits))}
            if (method=="jackknife") {
                resultx <- diversityresult(x, y, factor=factor1, level=levels[i], method="jackknife", index=index, digits=digits)
                result[i, 2] <- as.numeric(resultx$jack.estimate)
            }
        }
        if (sortit == T) {
            result2 <- result
            seq <- order(result[, 2])
            for (i in 1:m) {
                result[1:m, ] <- result2[seq, ]
            }
            rownames(result) <- rownames(result2)[seq]
        }
        return(result)

    }else{

        if (method == "jackknife") {stop("jackknife analysis problematic with two factors")}

        groups <- table(y[, factor1], y[, factor2])
        levels1 <- rownames(groups)
        levels2 <- colnames(groups)
        m1 <- length(levels1)
        m2 <- length(levels2)
        result <- array(NA, dim=c(m1, m2, 2))
        result[,,1] <- groups        
        dimnames(result) <- list(factor1=levels1, factor2=levels2, c("n", index))
        names(dimnames(result)) <- c(factor1, factor2, "")
        for (i in 1:m1) {
            for (j in 1:m2) {
                if (as.numeric(groups[i, j]) > 0) {
                    subs <- y[, factor1] == as.character(levels1[i])
                    x1 <- x[subs, , drop=F]
                    y1 <- y[subs, , drop=F]
                    result[i, j, 2] <- as.numeric(diversityresult(x1, y=y1, factor=factor2, level=levels2[j], method=method, index=index, digits=digits))
                }
            }
        }
        return(result)
    }
}

