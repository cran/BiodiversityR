`ensemble.simplified.categories` <- function(
    xcat=NULL, p=NULL, 
    filename=NULL, overwrite=TRUE, ...
)
{
    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(xcat,"RasterLayer") == F) {stop("parameter xcat is expected to be a RasterLayer")}
    if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
    if(is.null(filename) == T) {
        cat(paste("\n", "No new filename was provided", sep = ""))
        if (! require(tools)) {stop("tools package not available")}
        filename1 <- filename(xcat)
        extension1 <- paste(".", file_ext(filename1),sep="")
        extension2 <- paste("_new.", file_ext(filename1),sep="")
        filename <- gsub(pattern=extension1, replacement=extension2, x=filename1)
        cat(paste("\n", "New raster will be saved as:", sep = ""))
        cat(paste("\n", filename, "\n", sep=""))
    }
# get categories of presence points
    a <- randomPoints(xcat, n=10)
    TrainData <- prepareData(stack(xcat), p, b=a, factors=names(xcat), xy=FALSE)
    TrainData <- TrainData[TrainData[,"pb"]==1,]
    presence.categories <- levels(as.factor(TrainData[,names(xcat)]))
    presence.categories <- as.numeric(presence.categories)
    cat(paste("\n", "categories with presence points", "\n", sep = ""))
    print(presence.categories)

# get all categories of the layer
    all.categories <- freq(xcat)[,1]
    all.categories <- all.categories[is.na(all.categories) == F]

# categories without presence points
    new.categories <- all.categories[is.na(match(all.categories, presence.categories) )]
    cat(paste("\n", "categories without presence points", "\n", sep = ""))
    print(new.categories)

# outside category
    out.cat <- max(new.categories)
    cat(paste("\n", "categories without presence points will all be classified as: ", out.cat, "\n\n", sep = ""))
    replace.frame <- data.frame(id=new.categories, v=rep(out.cat, length(new.categories)))
    colnames(replace.frame)[2] <- names(xcat)
    new.x <- subs(xcat, replace.frame, by=1, which=2, subsWithNA=FALSE, filename=filename, overwrite=overwrite, ...)

    return(filename)
}
