`rankabuncomp` <-
function(x, y=NULL, factor=NULL,
    scale="abundance", scaledx=F, type="o", rainbow=T, legend=T, 
    xlim=c(1, max1), ylim=c(0, max2), 
    ...
) 
{
    groups <- table(y[,factor])
    levels <- names(groups)
    m <- length(groups)
    max1 <- max(diversitycomp(x, y, factor1=factor, index="richness", method="pooled")[,2])
    if (scaledx==T) {xlim<-c(0, 100)}
#    freq <- diversityresult(x, index="Berger", method="pooled")
    max2 <- max.2 <- 0
    for (i in 1:m) {
       if (scale=="abundance") {max.2 <- rankabundance(x, y, factor, levels[i])[1, "abundance"]}
       if (scale=="logabun") {max.2 <- rankabundance(x, y, factor, levels[i])[1, "abundance"]}
       if (scale=="proportion") {max.2 <- rankabundance(x, y, factor, levels[i])[1, "proportion"]}
       if (max.2 > max2) {max2 <- max.2}
    }
    if (scale=="accumfreq") {max2 <- 100}
    max2 <- as.numeric(max2)

    if (rainbow==F) {
        if (scale == "logabun" && all.equal(ylim, c(0, max2)) == T) {ylim <- c(1, max2)}
        rankabunplot(rankabundance(x, y, factor, levels[1]), scale=scale, scaledx=scaledx, type=type, labels=levels[1], xlim=xlim, ylim=ylim, pch=1, specnames=NULL, ...)
        for (i in 2:m) {
            rankabunplot(rankabundance(x, y, factor, levels[i]), addit=T, scale=scale, scaledx=scaledx, type=type, labels=levels[i], pch=i, specnames=NULL,...)
        }
        if (legend==T) {legend(graphics::locator(1), legend=levels, pch=c(1:m))}
    }else{
        grDevices::palette(rainbow(m))
        if (scale == "logabun" && all.equal(ylim, c(0, max2)) == T) {ylim <- c(1, max2)}
        rankabunplot(rankabundance(x, y, factor, levels[1]), scale=scale, scaledx=scaledx, type=type, labels=levels[1], xlim=xlim, ylim=ylim, col=1, pch=1, specnames=NULL,...)
        for (i in 2:m) {
            rankabunplot(rankabundance(x, y, factor, levels[i]), addit=T, scale=scale, scaledx=scaledx, type=type, labels=levels[i], col=i, pch=i, specnames=NULL,...)
        }
        if (legend==T) {legend(graphics::locator(1), legend=levels, pch=c(1:m), col=c(1:m))}
        grDevices::palette("default")
    }
}

