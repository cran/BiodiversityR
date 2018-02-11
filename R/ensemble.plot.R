`ensemble.plot` <- function(
    RASTER.species.name="Species001", RASTER.stack.name="base",
    plot.method=c("suitability", "presence", "count",
        "consensussuitability", "consensuspresence", "consensuscount", "consensussd"), 
    dev.new.width=7, dev.new.height=7,
    main=paste(RASTER.species.name, " ", plot.method, " for ", RASTER.stack.name, sep=""),
    positive.filters=c("grd"), negative.filters=c("xml"), 
    p=NULL, a=NULL,
    threshold=-1,
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    abs.breaks=6, abs.col=NULL,
    pres.breaks=6, pres.col=NULL,
    sd.breaks=9, sd.col=NULL,
    absencePresence.col=NULL, 
    count.col=NULL,
    maptools.boundaries=TRUE, maptools.col="dimgrey", ...
)
{

    plot.method <- match.arg(plot.method) 

    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (threshold < 0  && plot.method=="suitability") {
        if (is.null(p)==T || is.null(a)==T) {stop(paste("Please provide locations p and a to calculate thresholds", "\n", sep = ""))}
    }
    if(is.null(p) == F) {names(p) <- c("x", "y")}
    if(is.null(a) == F) {names(a) <- c("x", "y")}
#
# get raster files
# basic assumption is that different ensemble files are named as species_ENSEMBLE_1, species_ENSEMBLE_2, ... i.e. output from ensemble.batch

    species_focus <- RASTER.species.name
    if (gsub(".", "_", RASTER.species.name, fixed=T) != RASTER.species.name) {cat(paste("\n", "WARNING: species name (", RASTER.species.name, ") contains '.'", "\n\n", sep = ""))}
# 
   if (plot.method == "suitability") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//suitability", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to plot suitability as there are no files available for species: ", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    if (plot.method == "presence") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//presence", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to plot presence as there are no files available for species: ", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    if (plot.method == "count") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//count", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to plot counts as there are no files available for this species: ", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    if (plot.method == "consensussuitability") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//consensussuitability", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to plot consensus suitability as there are no files available for species: ", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    if (plot.method == "consensuspresence") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//consensuspresence", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to plot consensus presence as there are no files available for species: ", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    if (plot.method == "consensuscount") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//consensuscount", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to plot consensus counts as there are no files available for this species: ", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
    if (plot.method == "consensussd") {
        ensemble.files <- list.files(path=paste(getwd(), "//ensembles//consensussd", sep=""), pattern=species_focus, full.names=TRUE)
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to plot consensus standard deviations as there are no files available for this species: ", RASTER.species.name, "\n", sep = ""))
            return(NULL)
        }
    }
#
    RASTER.stack.name2 <- RASTER.stack.name
    if (gsub(".", "_", RASTER.stack.name, fixed=T) != RASTER.stack.name) {cat(paste("\n", "WARNING: title of stack (", RASTER.stack.name, ") contains '.'", "\n\n", sep = ""))}
    if (RASTER.stack.name != "") {
        ensemble.files <- ensemble.files[grepl(pattern=RASTER.stack.name, x=ensemble.files)]
        RASTER.stack.name2 <- paste("_", RASTER.stack.name, sep="")
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to plot as there are no raster files for this stack: ", RASTER.stack.name, "\n", sep = ""))
            return(NULL)
        }
    }
    for (i in 1:length(positive.filters)) {
        ensemble.files <- ensemble.files[grepl(pattern=positive.filters[i], x=ensemble.files)]
    }
    for (i in 1:length(negative.filters)) {
        ensemble.files <- ensemble.files[grepl(pattern=negative.filters[i], x=ensemble.files) == FALSE]
    }
    if (length(ensemble.files) < 1) {
        cat(paste("\n", "NOTE: not meaningful to plot as there are no raster files available for the specified filters", "\n", sep = ""))
        return(NULL)
    }
#
    cat(paste("\n", "Files used to create plots", "\n\n", sep = ""))
    print(ensemble.files)

    subtitle <- NULL
    threshold.mean <- threshold

    for (i in 1:length(ensemble.files)) {
        raster.focus <- raster::raster(ensemble.files[i])

        if (length(ensemble.files) > 1) {subtitle <- ensemble.files[i]}

        if (plot.method %in% c("suitability", "consensussuitability")) {
# thresholds apply to probabilities, also plot for probabilities
            raster.focus <- raster.focus / 1000
            raster.min <- raster::minValue(raster.focus)
            raster.max <- raster::maxValue(raster.focus)
            if (threshold.mean < 0) {
                eval1 <- NULL
                if (ncol(p) == 3) {p <- p[p[,1]==species_focus, c(2:3)]}
                if (ncol(a) == 3) {a <- a[a[,1]==species_focus, c(2:3)]}
                cat(paste("\n", "Evaluation of suitability raster layer at locations p and a", "\n", sep = ""))
                cat(paste("Note that threshold is only meaningful for calibration stack suitabilities", "\n\n", sep = ""))
                pres_consensus <- raster::extract(raster.focus, p)
                pres_consensus <- pres_consensus[is.na(pres_consensus)==F]
                abs_consensus <- raster::extract(raster.focus, a)
                abs_consensus <- abs_consensus[is.na(abs_consensus)==F]
                eval1 <- dismo::evaluate(p=pres_consensus, a=abs_consensus)
                print(eval1)
                threshold.mean <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres_consensus, Abs=abs_consensus)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(threshold.mean))
            }
            if (abs.breaks > 0) {
                seq1 <- round(seq(from=raster.min, to=threshold.mean, length.out=abs.breaks), 4)
                seq1 <- seq1[1:(abs.breaks-1)]
                seq1[-abs.breaks]
                seq1 <- unique(seq1)
                seq2 <- round(seq(from = threshold.mean, to = raster.max, length.out=pres.breaks), 4)
                seq2 <- unique(seq2)
                if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
                if (is.null(abs.col) == T) {abs.col <- grDevices::rainbow(n=length(seq1), start=0, end =1/6)}
                if (is.null(pres.col) == T) {pres.col <- grDevices::rainbow(n=length(seq2)-1, start=3/6, end=4/6)}
                raster::plot(raster.focus, breaks = c(seq1, seq2), col = c(abs.col, pres.col), colNA = NA, 
                    legend.shrink=0.8, cex.axis=0.8, main=main, sub=subtitle, ...)
            }else{
                seq1 <- NULL
                abs.col <- NULL
                seq2 <- round(seq(from = threshold.mean, to = raster.max, length.out=pres.breaks), 4)
                seq2 <- unique(seq2)
                if (is.null(pres.col) == T) {pres.col <- grDevices::rainbow(n=length(seq2)-1, start=3/6, end=4/6)}
                if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
                raster::plot(raster.focus, breaks = seq2, col = pres.col, colNA = NA, 
                    lab.breaks=seq2, legend.shrink=0.8, cex.axis=0.8, main=main, sub=subtitle, ...)
            }
        }

        if (plot.method %in% c("presence", "consensuspresence")) {
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            if (is.null(absencePresence.col) == T) {absencePresence.col <- c("grey", "green")}
            if (length(absencePresence.col) == 2) {
                raster::plot(raster.focus, breaks=c(0, 0.5, 1), col = absencePresence.col, colNA = NA, 
                    legend.shrink=0.6, cex.axis=0.8, lab.breaks=c("", "a", "p"), main=main, sub=subtitle, ...)
            }
            if (length(absencePresence.col) == 1) {
                raster::plot(raster.focus, breaks=c(0.5, 1), col = absencePresence.col, colNA = NA, 
                    legend.shrink=0.6, cex.axis=0.8, lab.breaks=c("a", "p"), main=main, sub=subtitle, ...)
            }
        }

        if (plot.method %in% c("count", "consensuscount")) {
            nmax <- raster::maxValue(raster.focus)
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            if (is.null(count.col) == T) {
                if (nmax > 3) {
                    count.col <- c("grey", "black", grDevices::rainbow(n=(nmax-2), start=0, end=1/3), "blue")
                }else{
                    count.col <- c("grey", grDevices::rainbow(n=nmax, start=0, end=1/3))
                }
            }
            seq1 <- seq(from=-1, to=nmax, by=1)
            if (length(count.col) == (length(seq1)-1)) {
                raster::plot(raster.focus, breaks=seq(from=-1, to=nmax, by=1), col=count.col, 
                    legend.shrink=0.8, cex.axis=0.8, main=main, sub=subtitle, ...)
            }else{
                raster::plot(raster.focus, breaks=c(0.1, seq(from=1, to=nmax, by=1)), col=count.col, 
                    lab.breaks=seq(from=0, to=nmax, by=1), legend.shrink=0.8, cex.axis=0.8, main=main, sub=subtitle, ...)
            }
        }

        if (plot.method == "consensussd") {
            sd.max <- raster::maxValue(raster.focus)
            seq1 <- seq(from = 0, to = sd.max, length.out = sd.breaks)
            if (is.null(sd.col) == T) {sd.col <- grDevices::rainbow(n = length(seq1)-1, start = 1/6, end = 4/6)}
            raster::plot(raster.focus, breaks=seq1, col=sd.col, 
                legend.shrink=0.8, cex.axis=0.8, main=main, sub=subtitle, ...)
       }

        if (maptools.boundaries == T) {
            data(wrld_simpl, package="maptools", envir=.BiodiversityR)
            maptools.wrld_simpl <- eval(as.name("wrld_simpl"), envir=.BiodiversityR)
            raster::plot(maptools.wrld_simpl, add=T, border=maptools.col)
        }

    }

    if (length(ensemble.files)==1 && plot.method %in% c("suitability", "consensussuitability")) {return(list(threshold=threshold.mean, breaks=c(seq1, seq2), 
        col=c(grDevices::rainbow(n=length(seq1), start = 0, end = 1/6), grDevices::rainbow(n=length(seq2)-1, start=3/6, end=4/6))))}

}
