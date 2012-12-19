`ensemble.weights` <- function(
    weights=c(0.9, 0.8, 0.7, 0.5), decay=1.5, multiply=TRUE, scale=TRUE, best=0
)
{
    weights <- as.array(weights)
    lw <- length(weights)
    lp <- sum(as.numeric(weights>0))
    if (best < 1) {best <- lp}
    if (lp < best) {
        cat(paste("\n","parameter best (",best,") larger than number of positive weights (", lp, ")","\n", sep=""))
        cat(paste("parameter best therefore ignored","\n", "\n", sep=""))
}
    original <- order(weights, decreasing=T)   
    decays <- as.numeric(rep(0,lw))
    decays[best] <- 1
    if (best > 1) {
        for (i in (best-1):1){
            decays[i] <- decays[i+1] * decay
        }
    }
#   check for equal ranks
    ranks <- rank(weights[original])
    ranklevels <- levels(as.factor(ranks))
    if (length(ranklevels) < lw) {
        cat(paste("\n","there were some ties in weights","\n", sep=""))
        cat(paste("average decay therefore used in some cases","\n", "\n", sep=""))        
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




