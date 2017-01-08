`ensemble.mean` <- function(
    RASTER.species.name="Species001", RASTER.stack.name="base",
    positive.filters=c("grd", "_ENSEMBLE_"), negative.filters=c("xml"), 
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
    abs.breaks=6, pres.breaks=6, sd.breaks=9,
    p=NULL, a=NULL,
    pt=NULL, at=NULL,
    threshold=-1,
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (threshold < 0) {
        if (is.null(p)==T || is.null(a)==T) {stop(paste("Please provide locations p and a to calculate thresholds", "\n", sep = ""))}
    }
    retest <- F
    if (is.null(pt)==F && is.null(at)==F) {
        if(identical(pt, p) == F || identical(at, a) == F)  {retest <- T}
    }
#
    if(is.null(p) == F) {names(p) <- c("x", "y")}
    if(is.null(a) == F) {names(a) <- c("x", "y")}
    if(is.null(pt) == F) {names(pt) <- c("x", "y")}
    if(is.null(at) == F) {names(at) <- c("x", "y")}
#
# avoid problems with non-existing directories
    dir.create("ensembles/consensussuitability", showWarnings = F)
    dir.create("ensembles/consensuscount", showWarnings = F)
    dir.create("ensembles/consensuspresence", showWarnings = F)
    dir.create("ensembles/consensussd", showWarnings = F)
    if(KML.out == T) {
        dir.create("kml", showWarnings = F)
        dir.create("kml/consensussuitability", showWarnings = F)
        dir.create("kml/consensuscount", showWarnings = F)
        dir.create("kml/consensuspresence", showWarnings = F)
        dir.create("kml/consensussd", showWarnings = F)
    }
#
# get ensemble input files
    species_focus <- RASTER.species.name
    if (gsub(".", "_", RASTER.species.name, fixed=T) != RASTER.species.name) {cat(paste("\n", "WARNING: species name (", RASTER.species.name, ") contains '.'", "\n\n", sep = ""))}
    ensemble.files <- list.files(path=paste(getwd(), "//ensembles//suitability", sep=""), pattern=species_focus, full.names=TRUE)
    if (length(ensemble.files) < 1) {
        cat(paste("\n", "NOTE: not meaningful to provide means as there are no raster files for this species:", RASTER.species.name, "\n", sep = ""))
        return(NULL)
    }
    RASTER.stack.name2 <- RASTER.stack.name
    if (gsub(".", "_", RASTER.stack.name, fixed=T) != RASTER.stack.name) {cat(paste("\n", "WARNING: title of stack (", RASTER.stack.name, ") contains '.'", "\n\n", sep = ""))}
    if (RASTER.stack.name != "") {
        ensemble.files <- ensemble.files[grepl(pattern=RASTER.stack.name, x=ensemble.files)]
        filename0 <- paste(species_focus, "_", RASTER.stack.name, sep="")
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to provide means as there are no raster files for this stack:", RASTER.stack.name, "\n", sep = ""))
            return(NULL)
        }
    }
    for (i in 1:length(positive.filters)) {
        ensemble.files <- ensemble.files[grepl(pattern=positive.filters[i], x=ensemble.files)]
    }
    for (i in 1:length(negative.filters)) {
        ensemble.files <- ensemble.files[grepl(pattern=negative.filters[i], x=ensemble.files) == FALSE]
    }
    if (length(ensemble.files) < 2) {
        cat(paste("\n", "NOTE: not meaningful to provide means as there are fewer than 2 ensemble files", "\n", sep = ""))
        return(NULL)
    }
    cat(paste("\n", "Files used to create mean ensemble", "\n\n", sep = ""))
    print(ensemble.files)
    ensemble.stack <- raster::stack(ensemble.files)
    cat(paste("\n", "RasterStack used to create mean ensemble", "\n\n", sep = ""))
    print(ensemble.stack)
#
    ensemble.mean <- raster::mean(ensemble.stack)
    ensemble.mean <- trunc(1.0 * ensemble.mean)
#    raster::setMinMax(ensemble.mean)
    names(ensemble.mean) <- filename0
    cat(paste("\n", "consensus mean suitability (truncated)", "\n\n", sep = ""))
    print(ensemble.mean)
    filename1 <- paste(getwd(), "//ensembles//consensussuitability//", filename0, sep="")
    raster::writeRaster(x=ensemble.mean, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
# avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble.mean, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- filename0
    raster::writeRaster(working.raster, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
# standard deviation
    ensemble.sd <- raster::calc(ensemble.stack, fun=sd)
    ensemble.sd <- trunc(ensemble.sd)
#    raster::setMinMax(ensemble.mean)
    names(ensemble.sd) <- filename0
    cat(paste("\n", "consensus standard deviation (truncated)", "\n\n", sep = ""))
    print(ensemble.sd)
    filename1 <- paste(getwd(), "//ensembles//consensussd//", filename0, sep="")
    raster::writeRaster(x=ensemble.sd, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble.sd, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- filename0
    raster::writeRaster(working.raster, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
#
    threshold.mean <- threshold
    if (threshold.mean < 0) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created mean ensemble raster layer at locations p and a", "\n", sep = ""))
        if (ncol(p) == 3) {p <- p[p[,1]==species_focus, c(2:3)]}
        if (ncol(a) == 3) {a <- a[a[,1]==species_focus, c(2:3)]}
        pres_consensus <- raster::extract(ensemble.mean, p)/1000
        pres_consensus <- pres_consensus[is.na(pres_consensus)==F]
        abs_consensus <- raster::extract(ensemble.mean, a)/1000
        abs_consensus <- abs_consensus[is.na(abs_consensus)==F]
        eval1 <- dismo::evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
        threshold.mean <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres_consensus, Abs=abs_consensus)
        cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
        print(as.numeric(threshold.mean))
        if(retest == T) {
            eval1 <- NULL
            cat(paste("\n", "Evaluation of created mean ensemble raster layer at locations pt and at", "\n\n", sep = ""))
            if (ncol(pt) == 3) {pt <- pt[pt[,1]==species_focus, c(2:3)]}
            if (ncol(at) == 3) {at <- at[at[,1]==species_focus, c(2:3)]}
            pres_consensus <- raster::extract(ensemble.mean, p)/1000
            pres_consensus <- pres_consensus[is.na(pres_consensus)==F]
            abs_consensus <- raster::extract(ensemble.mean, a)/1000
            abs_consensus <- abs_consensus[is.na(abs_consensus)==F]
            eval1 <- dismo::evaluate(p=pres_consensus, a=abs_consensus)
            print(eval1)
        }
    }

#
# 
    if (KML.out==T && raster::isLonLat(ensemble.mean)==F) {
        cat(paste("\n", "NOTE: not possible to generate KML files as Coordinate Reference System (CRS) is not longitude and latitude", "\n", sep = ""))
        KML.out <- FALSE
    }
    if (KML.out == T) {
        raster.min <- raster::minValue(ensemble.mean)
        raster.max <- raster::maxValue(ensemble.mean)
        seq1 <- round(seq(from=raster.min, to=threshold.mean, length.out=abs.breaks), 4)
        seq1 <- seq1[1:(abs.breaks-1)]
        seq1[-abs.breaks]
        seq1 <- unique(seq1)
        seq2 <- round(seq(from = threshold.mean, to = raster.max, length.out=pres.breaks), 4)
        seq2 <- unique(seq2)
        filename2 <- paste(getwd(), "//kml//consensussuitability//", filename0, sep="")
        raster::KML(ensemble.mean, filename=filename2, breaks = c(seq1, seq2), col = c(grDevices::rainbow(n=length(seq1), start=0, end =1/6), grDevices::rainbow(n=length(seq2)-1, start=3/6, end=4/6)), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE)
#
        sd.max <- raster::cellStats(ensemble.sd, stat='max')
        seq1 <- seq(from = 0, to = sd.max, length.out = sd.breaks)
        filename2b <- paste(getwd(), "//kml//consensussd//", filename0, sep="")
        raster::KML(ensemble.sd, filename=filename2b, col=grDevices::rainbow(n = length(seq1)-1, start = 1/6, end = 4/6), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE, breaks = seq1)
    }
#
# presence-absence maps based on the mean maps
    enspresence <- ensemble.mean >= 1000 * threshold.mean
    raster::setMinMax(enspresence)
    names(enspresence) <- filename0
    filename3 <- paste(getwd(), "//ensembles//consensuspresence//", filename0, sep="")
    raster::writeRaster(x=enspresence, filename=filename3, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(enspresence, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- filename0
    raster::writeRaster(working.raster, filename=filename3, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    if (KML.out == T) {
        filename4 <- paste(getwd(), "//kml//consensuspresence//", filename0, sep="")
        raster::KML(enspresence, filename=filename4, col=c("grey", "green"),
            colNA=0, blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE)
    }
#
# count maps: counting the number of ensembles predicting presence

    presence.files <- list.files(path=paste(getwd(), "//ensembles//presence", sep=""), pattern=species_focus, full.names=TRUE)
    if (RASTER.stack.name != "") {
        presence.files <- presence.files[grepl(pattern=RASTER.stack.name, x=presence.files)]
    }
    for (i in 1:length(positive.filters)) {
        presence.files <- presence.files[grepl(pattern=positive.filters[i], x=presence.files)]
    }
    for (i in 1:length(negative.filters)) {
        presence.files <- presence.files[grepl(pattern=negative.filters[i], x=presence.files) == FALSE]
    }

    ensemble.stack <- raster::stack(presence.files)
    cat(paste("\n", "RasterStack (presence-absence) used to create consensus ensemble (count of ensembles)", "\n\n", sep = ""))
    print(ensemble.stack)
    ensemble.count <- sum(ensemble.stack)
    raster::setMinMax(ensemble.count)
    names(ensemble.count) <- filename0
    filename5 <- paste(getwd(), "//ensembles//consensuscount//", filename0, sep="")
    raster::writeRaster(x=ensemble.count, filename=filename5, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble.count, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- filename0
    raster::writeRaster(working.raster, filename=filename5, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    if (KML.out == T) {
        filename6 <- paste(getwd(), "//kml//consensuscount//", filename0, sep="")
        nmax <- length(presence.files)
        if (nmax > 3) {
            raster::KML(ensemble.count, filename=filename6, col=c("grey", "black", grDevices::rainbow(n=(nmax-2), start=0, end=1/3), "blue"),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }else{
            raster::KML(ensemble.count, filename=filename6, col=c("grey", grDevices::rainbow(n=nmax, start=0, end=1/3)),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }
    }
    return(list(threshold=threshold.mean, call=match.call() ))
}

