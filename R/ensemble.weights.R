`ensemble.weights` <- function(
    weights=c(0.9, 0.8, 0.7, 0.5), 
    best=0, min.weight=0, 
    decay=1.0, interval.width=0.1,
    digits=4
)
{
    names.weights <- names(weights)
    weights <- as.numeric(weights)
    if(any(weights > 1.0)) {stop("Input weights are expected to be ranged between 0 and 1")}
    names(weights) <- names.weights
#   weights should not be negative
    if (min.weight < 0) {min.weight <- 0}
    weights[weights < min.weight] <- 0
    weights[is.na(weights)] <- 0
#
# special case if all weights are zero
    if (sum(weights) == 0) {return(weights)}
#
# select best weights
# ties are handled correctly
    lw <- length(weights)
    lp <- sum(as.numeric(weights > 0))
    if (best < 1) {best <- lp}
    if (lp < best) {best <- lp}
    weights.sorted <- sort(weights, decreasing=T)   
    min.best <- weights.sorted[best]
    weights[weights < min.best] <- 0
    weights.new <- weights
#
# create intervals
    intervals <- as.numeric(c(-1,min.weight))
    while ((max(intervals)+interval.width) <= 1.0) {intervals <- c(intervals, (max(intervals)+interval.width))}
#
    decays <- numeric(length=length(intervals))
    decays[1] <- 0
    decays[2] <- 1.0    
    if (decay < 1.0) {decay <- 1.0}
    if (length(decays) > 2) {
        for (i in 3: length(decays)) {decays[i] <- decays[i-1] * decay}
    }
#
    for (i in 1:length(weights)) {
        for (j in 1:(length(intervals)-1)) {
            if (intervals[j] <= weights[i] && weights[i] < intervals[j+1]) {weights.new[i] <- decays[j] * weights[i]}
        }
    }
    weights <- weights.new
# scaling to 1
    tot <- sum(weights)
    weights <- weights/tot
    weights <- round(weights, digits=digits)
#
    return(weights)
}
