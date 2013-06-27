`ensemble.weights` <- function(
    weights=c(0.9, 0.8, 0.7, 0.5), 
    min.weight=0,
    decay=1.5, multiply=TRUE, scale=TRUE, best=0
)
{
    names.weights <- names(weights)
    weights <- as.numeric(weights)
    names(weights) <- names.weights
#   weights should not be negative
    if (min.weight < 0) {min.weight <- 0}
    weights[weights < min.weight] <- 0
    lw <- length(weights)
    lp <- sum(as.numeric(weights > 0))
    if (best < 1) {best <- lp}
    if (lp < best) {best <- lp}
    original <- order(weights, decreasing=T)   
    decays <- as.numeric(rep(0,lw))
    decays[best] <- 1
    if (best > 1) {
        for (i in (best-1):1){
            decays[i] <- decays[i+1] * decay
        }
    }
#   check for equal ranks > 0
    ranks <- rank(weights[original])
    ranklevels <- levels(as.factor(ranks))
    if (lp < lw) {ranklevels <- ranklevels[ranklevels != as.character(min(as.numeric(ranklevels)))]}
    if (length(ranklevels) < lp  && decay != 1) {
# some ties in weights need to be replaced by mean decay
        for (i in 1:length(ranklevels)) {
            index <- ranks==ranklevels[i]
            meandecay <- mean(decays[index])
            decays[index] <- meandecay
        }
    }
    if (multiply==T) {
        for (i in 1:lw) {
            weights[original[i]] <- weights[original[i]] * decays [i]
        }
    }else{
        for (i in 1:lw) {
            weights[original[i]] <- decays [i]
        }        

    }
    if (scale==TRUE) {
        tot <- sum(weights)
        weights <- weights/tot
    }
    return(weights)
}



