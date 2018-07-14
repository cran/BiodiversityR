`ensemble.pairs` <- function(
    x=NULL, a=NULL, an=10000
)
{

# application of the graphics::pairs function including auxillary functions shown in the documentation (graphics 3.4.3)
# 
# raster package has a similar function 

## obtained from graphics::pairs example
## put histograms on the diagonal
    panel.hist <- function(x, ...){
        usr <- graphics::par("usr"); on.exit(graphics::par(usr))
        graphics::par(usr = c(usr[1:2], 0, 1.5) )
        h <- graphics::hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        graphics::rect(breaks[-nB], 0, breaks[-1], y, col = "green", ...)
    }

## obtained from graphics::pairs example
## modified to limit reduction of size to minimum 0.5 and show sign of correlation
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
        usr <- graphics::par("usr"); on.exit(graphics::par(usr))
        graphics::par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits = digits)[1]
## modified next 3 lines
        r1 <- max(0.5, abs(cor(x, y)))
        r2 <- cor(x, y)
        txt <- format(c(r2, 0.123456789), digits = digits)[1]
##
        txt <- paste0(prefix, txt)
        if(missing(cex.cor)) cex.cor <- 0.8/graphics::strwidth(txt)
## modified final line
##        text(0.5, 0.5, txt, cex = cex.cor * r)
        graphics::text(0.5, 0.5, txt, cex = cex.cor * r1)

    }

    if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
    if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
    if(is.null(a) == F) {names(a) <- c("x", "y")}
    if (is.null(a) == T) {a <- dismo::randomPoints(x[[1]], n=an, p=NULL, excludep=F)}
    x.data <- raster::extract(x, y=a)
    x <- as.data.frame(x.data)
    x <- na.omit(x)
    graphics::pairs(x, lower.panel=graphics::panel.smooth, diag.panel=panel.hist, upper.panel=panel.cor)
}

