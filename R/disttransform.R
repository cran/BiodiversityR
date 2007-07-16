`disttransform` <-
function(x,method="hellinger") {
    x <- as.matrix(x)
    METHODS <- c("hellinger","chord","profiles","chi.square","log","square","pa")
    method <- match.arg(method,METHODS)
    switch(method, hellinger = {
        x <- decostand(x,"hellinger")
    }, profiles = {
        x <- decostand(x,"total")
    }, chord = {
        x2 <- x^2
        rowtot <- apply(x2,1,sum)
        for (i in 1:length(rowtot)) {if (rowtot[i]==0) {rowtot[i] <- 1}}
        rowtot <- rowtot^0.5
        x <- x/rowtot
    }, chi.square = {
        x <- decostand(x,"chi.square")
    }, log = {
        x <- log(x+1)
    }, square = {
        x <- x^0.5
    }, pa = {
        x <- decostand(x,"pa")
    })
    x <- data.frame(x)
    return(x)
}

