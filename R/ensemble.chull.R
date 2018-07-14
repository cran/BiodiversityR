`ensemble.chull.create` <- function(
    x.pres=NULL, p=NULL, buffer.width=0.2,
    RASTER.format="raster", RASTER.datatype="INT1U", RASTER.NAflag=255,
    overwrite=TRUE,
    ...
)
{
#   if (! require(dismo)) {stop("Please install the dismo package")}
    if(is.null(x.pres) == T) {stop("value for argument x.pres is missing (RasterLayer object)")}
    if(inherits(x.pres, "RasterLayer") == F) {stop("x.pres is not a RasterLayer object")}
    x <- x.pres
    if (raster::maxValue(x) > 1) {
        cat(paste("Warning: base.raster has values larger than 1, hence does not provide presence-absence", sep=""))
    }
    if(is.null(p) == T) {stop("presence locations are missing")}
    names(p) <- c("x", "y")
#
# create convex hull around presence locations
# modified from red::map.habitat with option mcp=T
# modification includes creation of buffer around convex hull
#
#
    cat(paste("\n", "Creation of convex hull around presence locations", sep=""))
#
    vertices <- grDevices::chull(p)
    vertices <- c(vertices, vertices[1])
    vertices <- p[vertices, ]
    poly <- sp::Polygon(vertices)
    poly <- sp::Polygons(list(poly), 1)
    poly <- sp::SpatialPolygons(list(poly))
# modification for BiodiversityR
    maxdist1 <- max(raster::pointDistance(p, lonlat=F), na.rm=T)
    maxdist <- maxdist1 * buffer.width
    poly <- rgeos::gBuffer(poly, width=maxdist)
# modification ended
    patches <- raster::clump(x, gaps=F)
    selPatches <- raster::unique(raster::extract(patches, poly, df=T, weights=T)$clumps)
    selPatches <- selPatches[!is.na(selPatches)]
    allPatches <- raster::unique(patches)
    allPatches <- as.data.frame(cbind(allPatches, rep(0, length(allPatches))))
    names(allPatches) <- c("patches", "selected")
    allPatches[selPatches, 2] <- 1
    patches <- raster::subs(patches, allPatches)
    raster::setMinMax(patches)
#
    cat(paste("\n", "Buffer around convex hull of ", maxdist, " (", buffer.width, " * ", maxdist1, ", where ", maxdist1 , " is the maximum distance among presence locations)", "\n", sep=""))
#
# save
    raster.name <- names(x)
    names(patches) <- raster.name
    dir.create("ensembles/chull", showWarnings = F)
    filename1 <- paste("ensembles/chull/", raster.name, sep="")
    raster::writeRaster(patches, filename=filename1, progress='text', format=RASTER.format, overwrite=overwrite, datatype=RASTER.datatype, NAflag=RASTER.NAflag, ...)
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(patches, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- raster.name
    raster::writeRaster(working.raster, filename=filename1, progress='text', format=RASTER.format, overwrite=overwrite, datatype=RASTER.datatype, NAflag=RASTER.NAflag, ...)
#
    return(list(mask.layer=working.raster, convex.hull=poly))
}


`ensemble.chull.apply` <- function(
    x.spec=NULL, mask.layer=NULL, keep.old=T,
    RASTER.format="raster", RASTER.datatype="INT1U", RASTER.NAflag=255,
    overwrite=TRUE,
    ...
)
{
#   if (! require(dismo)) {stop("Please install the dismo package")}
    if(is.null(x.spec) == T) {stop("value for parameter x.spec is missing (RasterLayer object)")}
    if(inherits(x.spec, "RasterLayer") == F) {stop("x.spec is not a RasterLayer object")}
    x <- x.spec
    if(is.null(mask.layer) == T) {stop("value for parameter mask.layer is missing (RasterLayer object)")}
    if(inherits(mask.layer, "RasterLayer") == F) {stop("mask.layer is not a RasterLayer object")}
    if (raster::maxValue(mask.layer) > 1) {
        cat(paste("Warning: mask.layer has values larger than 1, hence does not provide presence mask layer", sep=""))
    }    
#    if (! require(tools)) {stop("tools package not available")}
    filename1 <- raster::filename(x)
#
    if (keep.old == T){
        extension1 <- paste(".", tools::file_ext(filename1), sep="")
        extension2 <- paste("_old.", tools::file_ext(filename1), sep="")
        filename2 <- gsub(pattern=extension1, replacement=extension2, x=filename1)
        raster::writeRaster(x, filename=filename2, overwrite=overwrite, ...)
    }
#
# cells that are 0 in the mask should become 0 in the output file
    masked.x <- raster::mask(x, mask=mask.layer, maskvalue=0, updatevalue=0)
#
    raster.name <- names(x)
    names(masked.x) <- raster.name
    raster::writeRaster(masked.x, filename=filename1, progress='text', format=RASTER.format, overwrite=overwrite, datatype=RASTER.datatype, NAflag=RASTER.NAflag, ...)
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(masked.x, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- raster.name
    raster::writeRaster(working.raster, filename=filename1, progress='text', format=RASTER.format, overwrite=overwrite, datatype=RASTER.datatype, NAflag=RASTER.NAflag, ...)
#
    if (keep.old == T) {
        old.raster <- raster::raster(filename2)
        return(list(masked.raster=masked.x, old.raster=old.raster))
    }else{
        return(masked.x)
    }
}

