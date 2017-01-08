`renyiaccumresult` <- function(
    x, y=NULL, factor, level,
    scales=c(0, 0.25, 0.5, 1, 2, 4, 8, Inf),
    permutations=100,...) 
{
    if(is.null(y) == F) {

        if((factor %in% names(y)) == F) {stop("specified factor '", factor, "' is not a variable of the environmental data frame")}
        if(is.factor(y[, factor]) == F) {stop("specified factor '", factor, "' is not a factor")} 
        levels1 <- as.character(levels(as.factor(as.character(y[, factor]))))
        if((level %in% levels1) == F) {stop("specified level '", level, "' is not an available factor level")}

        subs <- y[, factor] == level
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q] <- F}
         }

         x <- x[subs,,drop=F]
         freq <- apply(x,2,sum)
         subs <- freq>0
         x <- x[, subs, drop=F]
    }

    result <- renyiaccum(x, scales=scales, permutations=permutations, ...)
    return(result)
}


