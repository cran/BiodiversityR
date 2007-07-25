`residualssurface` <-
function(model,data,x,y,gam=F,npol=2,plotit=T,filled=F,bubble=F) {
    if (!require(spatial)) {stop("Requires package spatial")}
    if (!require(akima)) {stop("Requires package akima")}
    res <- na.omit(residuals(model))
    xpos <- data[,x]
    ypos <- data[,y]
    if (gam==T) {
        if (!require(mgcv)) {stop("Requires package mgcv")}
        result <- gam(res~s(xpos)+s(ypos),family=gaussian)
    }else{
        result <- surf.ls(npol,xpos,ypos,res)
    }
    if (plotit==T) {
        fitted <- fitted(result)
        if (filled==F) {
            contour(interp(xpos,ypos,fitted,duplicate="mean"),lwd=2)
        }else{
            filled.contour(interp(xpos,ypos,fitted,duplicate="mean"),color=terrain.colors)
        }
        if (bubble==T && filled==F) {
            res2 <- abs(res)
            for (i in 1:length(res)) {
                if (res[i] < 0) {
                    res[i] <- NA
                }else{
                    res2[i] <- NA
                }
            }
            symbols(xpos,ypos,circles=res,add=T)
            symbols(xpos,ypos,squares=res2,add=T)
        }
        if (bubble==F && filled==F) {
            points(xpos,ypos)
        }
    }
    return(result)
}

