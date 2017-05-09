`ordisymbol` <-
function(ordiplot, y, factor, col=1, colors=TRUE, pchs=TRUE,
    rainbow=TRUE, heat.colors=FALSE, terrain.colors=FALSE, 
    topo.colors=FALSE, cm.colors=FALSE, 
    legend=TRUE, legend.x="topleft", legend.ncol=1, ...) 
{
    ordiscores <- scores(ordiplot, display="sites")
    groups <- table(y[,factor])
    m <- length(groups)
    if (m > 25) {
        warning("Symbol size was kept constant as there were more than 25 categories (> number of symbols that are currently used in R)")
        colors <- TRUE
        pchs <- FALSE
    }
    levels <- names(groups)
    if (rainbow == T) {
        grDevices::palette(rainbow(m))
        colors <- TRUE
    }
    if (heat.colors == T) {
        grDevices::palette(heat.colors(m))
        colors <- TRUE
    }
    if (terrain.colors == T) {
        grDevices::palette(terrain.colors(m))
        colors <- TRUE
    }
    if (topo.colors == T) {
        grDevices::palette(topo.colors(m))
        colors <- TRUE
    }
    if (cm.colors == T) {
        grDevices::palette(topo.colors(m))
        colors <- TRUE
    }
    for (i in 1:m) {
        subs <- y[,factor]==levels[i]
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q]<-F}
        }
        scores <- ordiscores[subs,,drop=F]
        if (colors==T && pchs==T) {
            graphics::points(scores[,1], scores[,2], pch=i, col=i,...)
        }
        if (colors==T && pchs==F) {
            graphics::points(scores[,1], scores[,2], pch=19, col=i,...)
        }
        if (colors == F) {
            graphics::points(scores[,1], scores[,2], pch=i, col=col,...)
        }
    }
    if (legend==T) {
        if (colors==T && pchs==T) {legend(x=legend.x, legend=levels, pch=c(1:m), col=c(1:m), ncol=legend.ncol)}
        if (colors==T && pchs==F) {legend(x=legend.x, legend=levels, pch=rep(19, m), col=c(1:m), ncol=legend.ncol)}
        if (colors == F) {legend(x=legend.x, legend=levels, pch=c(1:m))}
    }
    grDevices::palette("default")
}

