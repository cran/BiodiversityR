`ordispider2` <-
function (ordiplot, groups, display = "sites", w = weights(ordiplot, display), 
    show.groups, return.index=F,...) 
{
    if (inherits(ordiplot, "cca") && missing(groups)) {
        lc <- scores(ordiplot, display = "lc", ...)
        wa <- scores(ordiplot, display = "wa", ...)
        segments(lc[, 1], lc[, 2], wa[, 1], wa[, 2], ...)
        return(invisible())
    }
    ordispidercentres <- function (ordiplot, groups, display = "sites", w = weights(ordiplot, display), 
        show.groups) 
    {
        pts <- scores(ordiplot, display = display)
        w <- eval(weights(ordiplot, display))
        if (length(w) == 1) 
            w <- rep(1, nrow(pts))
        if (is.null(w)) 
            w <- rep(1, nrow(pts))
        if (!missing(show.groups)) {
            take <- groups %in% show.groups
            pts <- pts[take, , drop = FALSE]
            groups <- groups[take]
            w <- w[take]
        }
        out <- seq(along = groups)
        inds <- names(table(groups))
        result <- array(NA,dim=c(length(inds),2))
        rownames(result) <- inds
        for (is in inds) {
            gr <- out[groups == is]
            if (length(gr) > 1) {
                X <- pts[gr, ]
                W <- w[gr]
                ave <- apply(X, 2, weighted.mean, w = W)
                result[is,1] <- ave[1]
                result[is,2] <- ave[2]
            }
        }
        return(result)
    }
    pts <- scores(ordiplot, display = display, ...)
    w <- eval(w)
    if (length(w) == 1) 
        w <- rep(1, nrow(pts))
    if (is.null(w)) 
        w <- rep(1, nrow(pts))
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
        w <- w[take]
    }
    out <- seq(along = groups)
    inds <- names(table(groups))
    centres <- ordispidercentres(ordiplot=ordiplot,groups=groups)
    index <- rep(FALSE,length(out))
    g <- length(inds)
    for (i in 1:length(index)) {   
        mat1 <- rbind(centres,pts[i,])
        distances <- vegdist(mat1,"euc")
        distances <- as.matrix(distances)
        distances <- distances[1:g,g+1]
        ind2 <- groups[i]==names(distances)
        if(distances[ind2]==min(distances)) {index[i] <- TRUE}
    }
    groups <- groups[index]
    pts <- pts[index,]
    w <- w[index]
    out <- seq(along = groups)
    inds <- names(table(groups))
    for (is in inds) {
        gr <- out[groups == is]
        if (length(gr) > 1) {
            X <- pts[gr, ]
            W <- w[gr]
            ave <- apply(X, 2, weighted.mean, w = W)
            segments(ave[1], ave[2], X[, 1], X[, 2], ...)
        }
    }
    if (return.index ==T){
        cat("Points that are closest to their own ordispider centre\n")
        return(index)
    }
    invisible()
}



