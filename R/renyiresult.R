`renyiresult` <- function(
    x, y=NULL, factor, level, method="all",
    scales=c(0, 0.25, 0.5, 1, 2, 4, 8, Inf), evenness=FALSE, ...) 
{
    if (is.null(y) == F) {
        if((factor %in% names(y)) == F) {stop("specified factor '", factor, "' is not a variable of the environmental data frame")}
        if(is.factor(y[, factor]) == F) {stop("specified factor '", factor, "' is not a factor")} 
        levels1 <- as.character(levels(as.factor(as.character(y[, factor]))))
        if((level %in% levels1) == F) {stop("specified level '", level, "' is not an available factor level")}

        subs <- y[,factor] == as.character(level)
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q]<-F}
        }
        x <- x[subs,, drop=F]
        freq <- apply(x, 2, sum)
        subs <- freq > 0
        x <- x[, subs, drop=F]
    }
    if(method=="all") {x <- t(as.matrix(apply(x,2,sum)))}
    result <- renyi(x, scales=scales,...)
    if (attributes(result)$class[2] == "numeric") {
        result <- data.frame(t(as.matrix(result)))
        rownames(result) <- "all"
        colnames(result) <- scales        
    }
    if (evenness == T) {result[,] <- result[,]-renyi(x,scales=c(0))}
    return(result)
}



