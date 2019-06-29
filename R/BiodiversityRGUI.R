`BiodiversityRGUI` <-
function(changeLog=FALSE, backward.compatibility.messages=FALSE)
{

    if (backward.compatibility.messages == T) {
        cat(paste("\n", "Notes on backward compatiblity from BiodiversityR version 2.8.0", "\n"))
        cat(paste("\n", "In prior versions, function ensemble.calibrate.models was function ensemble.test"))
        cat(paste("\n", "For possible backward compatibility assign ensemble.test <- ensemble.calibrate.models"))
        cat(paste("\n", "In prior versions, function ensemble.calibrate.weights was function ensemble.test.splits"))
        cat(paste("\n", "In prior versions, slot ensemble.calibrate.weights$AUC.table was ensemble.calibrate.weights$table"))
        cat(paste("\n", "In prior versions, argument SSB.reduce was CIRCLES.at"))

        cat(paste("\n\n", "(The earlier name of ensemble.test originated from the first [2012] version of ensemble suitability"))
        cat(paste("\n", "modelling where both ensemble.raster and ensemble.test internally calibrated and evaluated [tested]"))
        cat(paste("\n", "models, but only ensemble.raster went ahead with creating suitability raster layers.)", "\n\n\n"))
    }

    if (changeLog == T) {BiodiversityR.changeLog()}

    if (! requireNamespace("vegan")) {stop("Please install the vegan package")}
    if (! requireNamespace("rgl")) {stop("Please install the rgl package")}    
    if (! requireNamespace("vegan3d")) {stop("Please install the vegan3d package")}
    if (! requireNamespace("dismo")) {stop("Please install the dismo package")}

    options(Rcmdr=list(etc=file.path(path.package(package="BiodiversityR"),
        "etc"), sort.names=FALSE))

    if ("Rcmdr" %in% .packages()) {
        stop("R commander should not have been loaded yet")
    }else{
        Rcmdr::Commander()
    }
}

