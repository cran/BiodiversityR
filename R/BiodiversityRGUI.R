`BiodiversityRGUI` <-
function()
{
#    if (! require(vegan)) {stop("Please install the vegan package")}
    if (! requireNamespace("vegan")) {stop("Please install the vegan package")}
    if (! requireNamespace("dismo")) {stop("Please install the dismo package")}
    options(Rcmdr=list(etc=file.path(path.package(package="BiodiversityR"),
        "etc"), sort.names=FALSE))
    if ("Rcmdr" %in% .packages()) {
        stop("R commander should not have been loaded yet")
    }else{
        Rcmdr::Commander()
    }
}

