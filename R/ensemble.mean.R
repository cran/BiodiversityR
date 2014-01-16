`ensemble.mean` <- function(
    RASTER.species.name="Species001", RASTER.stack.name="base",
    RASTER.extension="grd", 
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
    p=NULL, a=NULL,
    pt=NULL, at=NULL,
    threshold.method="spec_sens", threshold.sensitivity=0.9
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(p)==T || is.null(a)==T) {stop(paste("Please provide locations p and a to calculate thresholds", "\n", sep = ""))}
    retest <- F
    if (is.null(pt)==F && is.null(at)==F) {
        if(identical(pt, p) == F || identical(at, a) == F)  {retest <- T}
    }
#
# avoid problems with non-existing directories
    dir.create("ensembles/count", showWarnings = F)
    dir.create("ensembles/presence", showWarnings = F)
    if(KML.out == T) {
        dir.create("kml", showWarnings = F)
        dir.create("kml/count", showWarnings = F)
        dir.create("kml/presence", showWarnings = F)
    }


# get ensemble files
# basic assumption is that different ensemble files are named as species_ENSEMBLE_1, species_ENSEMBLE_2, ... i.e. output from ensemble.batch

    species_focus <- RASTER.species.name
    ensemble.files <- list.files(path=paste(getwd(), "//ensembles", sep=""), pattern=paste(species_focus, "_ENSEMBLE_", sep=""), full.names=TRUE)
    ensemble.files <- ensemble.files[grepl(pattern=RASTER.stack.name, x=ensemble.files)]
    ensemble.files <- ensemble.files[grepl(pattern=RASTER.extension, x=ensemble.files)]

    if (length(ensemble.files) < 2) {
        cat(paste("\n", "NOTE: not meaningful to provide means as there are fewer than 2 ensemble files", "\n", sep = ""))
        return(NULL)
    }
    ensemble.stack <- stack(ensemble.files)
    cat(paste("\n", "RasterStack used to create mean ensemble", "\n\n", sep = ""))
    print(ensemble.stack)
    ensemble.mean <- mean(ensemble.stack)
    if (is.na(ensemble.mean@crs) == T) {
        cat(paste("\n", "Problem with CRS of mean ensemble", "\n", sep = ""))
        return()
    }
    ensemble.mean <- trunc(ensemble.mean)
    setMinMax(ensemble.mean)
    names(ensemble.mean) <- paste(species_focus, "_MEAN", sep="")
    filename1 <- paste(getwd(), "//ensembles//", species_focus, "_MEAN_", RASTER.stack.name, sep="")
    writeRaster(x=ensemble.mean, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
    eval1 <- NULL
    cat(paste("\n", "Evaluation of created mean ensemble raster layer at locations p and a", "\n", sep = ""))
    pres_consensus <- extract(ensemble.mean, p)
    abs_consensus <- extract(ensemble.mean, a)
    eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
    print(eval1)
    threshold.mean <- threshold(eval1, sensitivity=threshold.sensitivity)[[threshold.method]]
    if(retest == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created mean ensemble raster layer at locations pt and at", "\n\n", sep = ""))
        pres_consensus <- extract(ensemble.mean, pt)
        abs_consensus <- extract(ensemble.mean, at)
        eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
    }
#
    if (KML.out == T) {
        seq1 <- seq(from = 0, to = threshold.mean, length.out = 10)
        seq2 <- seq(from = threshold.mean, to = 1000, length.out = 11)
        filename2 <- paste(getwd(), "//KML//", species_focus, "_MEAN_", RASTER.stack.name, sep="")
        KML(ensemble.mean, filename=filename2, col = c(rainbow(n = 10, start = 0, end = 1/6), rainbow(n = 10, start = 3/6, end = 4/6)), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE, breaks = c(seq1, seq2))
    }
#
# presence-absence maps based on the mean maps
    enspresence <- ensemble.mean > threshold.mean
    if (is.na(enspresence@crs) == T) {enspresence@crs <- ensemble.mean@crs}
    enspresence <- trunc(enspresence)
    setMinMax(enspresence)
    names(enspresence) <- paste(species_focus, "_MEAN_presence", sep="")
    filename3 <- paste(getwd(), "//ensembles//presence//", species_focus, "_MEAN_", RASTER.stack.name, sep="")
    writeRaster(x=enspresence, filename=filename3, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if (KML.out == T) {
        filename4 <- paste(getwd(), "//kml//presence//", species_focus, "_MEAN_", RASTER.stack.name, sep="")
        KML(enspresence, filename=filename4, col=c("grey", "green"),
            colNA=0, blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE)
    }
#
# count maps: counting the number of ensembles predicting presence
    presence.files <- list.files(path=paste(getwd(), "//ensembles//presence", sep=""), pattern=paste(species_focus, "_ENSEMBLE_", sep=""), full.names=TRUE)
    presence.files <- presence.files[grepl(pattern=RASTER.stack.name, x=presence.files)]
    presence.files <- presence.files[grepl(pattern=RASTER.extension, x=presence.files)]

    ensemble.stack <- stack(presence.files)
    cat(paste("\n", "RasterStack (presence-absence) used to create mean ensemble (count)", "\n\n", sep = ""))
    print(ensemble.stack)
    ensemble.count <- sum(ensemble.stack)
    if (is.na(ensemble.count@crs) == T) {ensemble.count@crs <- ensemble.mean@crs}
    ensemble.count <- trunc(ensemble.count)
    setMinMax(ensemble.count)
    names(ensemble.count) <- paste(species_focus, "_MEAN_count", sep="")
    filename5 <- paste(getwd(), "//ensembles//count//", species_focus, "_MEAN_", RASTER.stack.name, sep="")
    writeRaster(x=ensemble.count, filename=filename5, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if (KML.out == T) {
        filename6 <- paste(getwd(), "//kml//count//", species_focus, "_MEAN_", RASTER.stack.name, sep="")
        nmax <- length(presence.files)
        if (nmax > 3) {
            KML(ensemble.count, filename=filename6, col=c("grey", rainbow(n=(nmax-1), start=0, end=1/3), "blue"),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }else{
            KML(ensemble.count, filename=filename6, col=c("grey", rainbow(n=nmax, start=0, end=1/3)),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }
    }
}

