`ensemble.chull.create` <- function(
    x.pres=NULL, p=NULL, buffer.width=0.2, buffer.maxmins=FALSE, lonlat.dist=FALSE,
    poly.only=FALSE,
    RASTER.format="GTiff", RASTER.datatype="INT1U", RASTER.NAflag=255,
    overwrite=TRUE,
    ...
)
{
#   if (! require(dismo)) {stop("Please install the dismo package")}
    if (poly.only == F) {
        if(is.null(x.pres) == T) {stop("value for argument x.pres is missing (RasterLayer object)")}
        if(inherits(x.pres, "RasterLayer") == F) {stop("x.pres is not a RasterLayer object")}
        x <- x.pres
        if (raster::maxValue(x) > 1) {
            cat(paste("Warning: base.raster has values larger than 1, hence does not provide presence-absence", sep=""))
        }
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
    if (buffer.maxmins == FALSE) {
        maxdist1 <- max(raster::pointDistance(p, lonlat=F), na.rm=T)
        maxdist <- maxdist1 * buffer.width
        if (lonlat.dist == TRUE) {maxdist2 <- max(raster::pointDistance(p, lonlat=T), na.rm=T)}
    }else{
        point.dists <- raster::pointDistance(p, lonlat=F, allpairs=F)
        point.dists[point.dists == 0] <- NA
        maxdist1 <- max(apply(point.dists, 2, min, na.rm=T))
        if (lonlat.dist == TRUE) {
            pres.distances <- array(dim=c(nrow(p), nrow(p)))
            for (i in 1:nrow(p)) {
                pres.distances[, i] <- geosphere::distGeo(p[], p[i, ])
            }
            pres.distances <- data.frame(pres.distances)
            pres.distances[pres.distances == 0] <- NA
            min.distances <- apply(pres.distances, 2, min, na.rm=T)
            maxdist2 <- max(min.distances)
        }
        maxdist <- maxdist1 * buffer.width
    }
# changed Dec 2022 to remove dependency on rgeos
#    poly <- rgeos::gBuffer(poly, width=maxdist)
    poly.sf <- sf::st_as_sf(poly)
    poly <- sf::st_buffer(poly.sf, dist=maxdist)
#
    if (buffer.maxmins == FALSE) {
        cat(paste("\n", "Buffer around convex hull of ", maxdist, " (", buffer.width, " * ", maxdist1, ", where ", maxdist1 , " is the maximum distance among presence locations)", "\n", sep=""))
        if (lonlat.dist == TRUE) {cat(paste("This maximum distance corresponds to a distance in km of: ", maxdist2/1000, "\n"))}
    }else{
        cat(paste("\n", "Buffer around convex hull of ", maxdist, " (", buffer.width, " * ", maxdist1, ", where ", maxdist1 , " is the maximum of the distances to the closest neighbour for each presence location)", "\n", sep=""))
        if (lonlat.dist == TRUE) {cat(paste("This maximum distance corresponds to a distance in km of: ", maxdist2/1000, "\n"))}
    }
    
    if (poly.only == TRUE) {return(list(convex.hull=poly))}    
    
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
# save
    raster.name <- names(x)
    names(patches) <- raster.name
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/chull", showWarnings = F)
# modified to save the patches as 'mask'...
    filename1 <- paste("ensembles/chull/", "mask", sep="")
    raster::writeRaster(patches, filename=filename1, progress='text', format=RASTER.format, overwrite=overwrite, datatype=RASTER.datatype, NAflag=RASTER.NAflag, ...)
#  avoid possible problems with saving of names of the raster layers
# no longer used with default format of GTiff since DEC-2022
#    raster::writeRaster(patches, filename="working.grd", overwrite=T)
#    working.raster <- raster::raster("working.grd")
#    names(working.raster) <- raster.name
#    raster::writeRaster(working.raster, filename=filename1, progress='text', format=RASTER.format, overwrite=overwrite, datatype=RASTER.datatype, NAflag=RASTER.NAflag, ...)
#
    return(list(mask.layer=patches, convex.hull=poly))
}


`ensemble.chull.apply` <- function(
    x.spec=NULL, mask.layer=NULL, keep.old=T,
    RASTER.format="GTiff", RASTER.datatype="INT1U", RASTER.NAflag=255,
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
# modified in DEC-2022 to create a new directory rather than overwriting the presence file
#    if (keep.old == T){
#        extension1 <- paste(".", tools::file_ext(filename1), sep="")
#        extension2 <- paste("_old.", tools::file_ext(filename1), sep="")
#        filename2 <- gsub(pattern=extension1, replacement=extension2, x=filename1)
#        raster::writeRaster(x, filename=filename2, overwrite=overwrite, ...)
#    }
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/chull", showWarnings = F)
    filename3 <- paste0(getwd(), "/ensembles/chull/", basename(filename1))
#
# cells that are 0 in the mask should become 0 in the output file
    masked.x <- raster::mask(x, mask=mask.layer, maskvalue=0, updatevalue=0)
#
    raster.name <- names(x)
    names(masked.x) <- raster.name
    raster::writeRaster(masked.x, filename=filename3, progress='text', format=RASTER.format, overwrite=overwrite, datatype=RASTER.datatype, NAflag=RASTER.NAflag, ...)
#  avoid possible problems with saving of names of the raster layers
# no longer used with default format of GTiff since DEC-2022
#    raster::writeRaster(masked.x, filename="working.grd", overwrite=T)
#    working.raster <- raster::raster("working.grd")
#    names(working.raster) <- raster.name
#    raster::writeRaster(working.raster, filename=filename1, progress='text', format=RASTER.format, overwrite=overwrite, datatype=RASTER.datatype, NAflag=RASTER.NAflag, ...)
#
    if (keep.old == T) {
#        old.raster <- raster::raster(filename2)
        return(list(masked.raster=masked.x, old.raster=x))
    }else{
        return(masked.x)
    }
}

`ensemble.chull.buffer.distances` <- function(
    p=NULL, buffer.maxmins=FALSE, lonlat.dist=FALSE
)
{
    names(p) <- c("x", "y")

    if (buffer.maxmins==FALSE) {

# maximum
        point.dists <- raster::pointDistance(p, lonlat=F)
        max.distances <- apply(point.dists, 2, max, na.rm=T)
        max.1a <- which.max(max.distances)
        max.1b <- which.max(point.dists[max.1a, ])
        max.1c <- max(max.distances)
    
        cat(paste("Maximum distance is between locations ", max.1a, " and ", max.1b, " with distance in native coordinates of ",  max.1c, "\n", sep=""))

        if (lonlat.dist == TRUE) {
            pres.distances <- array(dim=c(nrow(p), nrow(p)))
            for (i in 1:nrow(p)) {
                pres.distances[, i] <- geosphere::distGeo(p[], p[i, ])
            }
            pres.distances <- data.frame(pres.distances)
            pres.distances[pres.distances == 0] <- NA
            max.distances <- apply(pres.distances, 2, max, na.rm=T)
            max.1a <- which.max(max.distances)
            max.1b <- which.max(pres.distances[max.1a, ])
            max.1c <- max(max.distances)/1000   
            cat(paste("Maximum distance is between locations ", max.1a, " and ", max.1b, " with distance in km of ",  max.1c, "\n", sep=""))
        }

        result <- c(max.1a, max.1b, max.1c)
        names(result) <- c("location.1", "location.2", "distance")

    }else{
    
        point.dists <- raster::pointDistance(p, lonlat=F)
        point.dists[point.dists == 0] <- NA
        min.distances <- apply(point.dists, 2, min, na.rm=T)
        max.1a <- which.max(min.distances)
        min.1b <- which.min(point.dists[max.1a, ])
        max.1c <- max(min.distances)
    
        cat(paste("Maximum closest neighbour distance is between locations ", max.1a, " and ", min.1b, " with distance in native coordinates of ",  max.1c, "\n", sep=""))

        if (lonlat.dist == TRUE) {
            pres.distances <- array(dim=c(nrow(p), nrow(p)))
            for (i in 1:nrow(p)) {
                pres.distances[, i] <- geosphere::distGeo(p[], p[i, ])
            }
            pres.distances <- data.frame(pres.distances)
            pres.distances[pres.distances == 0] <- NA
            min.distances <- apply(pres.distances, 2, min, na.rm=T)
            max.1a <- which.max(min.distances)
            min.1b <- which.min(pres.distances[max.1a, ])
            max.1c <- max(min.distances)/1000   
            cat(paste("Maximum closest neighbour distance is between locations ", max.1a, " and ", min.1b, " with distance in km of ",  max.1c, "\n", sep=""))      
        }    

        result <- c(max.1a, min.1b, max.1c)
        names(result) <- c("location.1", "location.2", "distance")
    
    }    

    return(result)
    
}    


`ensemble.chull.MSDM` <- function(
    p=NULL, a=NULL, species.name=NULL,
    suit.file=NULL, suit.divide=1000, MSDM.dir = NULL,
    method = "BMCP", threshold = "spec_sens", 
    buffer = "species_specific"
)
{
    if(is.null(p) == T) {stop("Provide p (presence observations)")}
    if(ncol(p) != 2) {stop("p (presence observations) is expected to have 2 columns (x and y coordinates)")}
    if(is.null(a) == T) {stop("Provide a (absence or background observations)")}
    if(ncol(a) != 2) {stop("a (absence observations) is expected to have 2 columns (x and y coordinates)")}
    if(is.null(species.name) == T) {stop("Provide species name")}
    if(file.exists(suit.file) == F) {stop("Suitability file does not exist")}  
    if(dir.exists(MSDM.dir) == F) {stop("MSDM directory does not exist")}  

    name.OK <- gsub(species.name, pattern=" ", replacement="_")
    p <- data.frame(sp=rep(name.OK, nrow(p)), p)
    a <- data.frame(sp=rep(name.OK, nrow(a)), a)
    names(p) <- names(a) <- c("sp", "x", "y")
    suit.raster <- raster::raster(suit.file)
    prob.raster <- suit.raster/1000
    names(prob.raster) <- name.OK
    prob.file <- paste(MSDM.dir, "/", name.OK, ".tif", sep="")
    raster::writeRaster(prob.raster, filename=prob.file, overwrite=TRUE)
    cat(paste("created file: ", prob.file, "\n", sep=""))

    MSDM.script <- paste("MSDM_Posteriori(records=M.out$records, absences=M.out$absences, ", "x='x', y='y', sp='sp',", "dirraster='", MSDM.dir, "', dirsave='", MSDM.dir, "', ", "method='", method, "', buffer='", buffer, "')", sep="")

    line1 <- paste("MSDM_Posteriori(records=M.out$records, absences=M.out$absences,", sep="")
    line2 <- paste("x='x', y='y', sp='sp',", sep="")
    line3 <- paste("dirraster='", MSDM.dir, "', dirsave='", MSDM.dir, "',", sep="") 
    line4 <- paste("method='", method, "', buffer='", buffer, "')", sep="")

    return(list(records=p, absences=a, script=MSDM.script,
        line1=line1, line2=line2, line3=line3, line4=line4))
}
