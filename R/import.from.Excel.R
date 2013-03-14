if(.Platform$OS.type == "windows") {

`import.from.Excel` <-
function(file=file.choose(), data.type="community", sheet=NULL, sitenames="sites", column="species", value="abundance", factor="", level="", cepnames=FALSE) {
    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- odbcConnectExcel(file)
    if (is.null(data.type) == TRUE) {data.type <- sheet}
    TYPES <- c("community", "environmental", "stacked")
    data.type <- match.arg(data.type, TYPES)
    if (is.null(sheet) == TRUE) {sheet <- data.type}
    if (data.type == "stacked") {
        stackeddata <- sqlFetch(dataplace,sheet)  
        result <- makecommunitydataset(stackeddata,row=sitenames,column=column,value=value,factor=factor,level=level)
    }else{
         result <- sqlFetch(dataplace,sheet,rownames=sitenames)
    }
    close(dataplace)
    rownames(result) <- make.names(rownames(result),unique=T)
    if (cepnames == TRUE) {
        colnames(result) <- make.cepnames(colnames(result))
    }else{
        colnames(result) <- make.names(colnames(result),unique=T)
    }
    return(result)
}

`import.from.Excel2007` <-
function(file=file.choose(), data.type="community", sheet=NULL, sitenames="sites", column="species", value="abundance", factor="", level="", cepnames=FALSE) {
    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- odbcConnectExcel2007(file)
    if (is.null(data.type) == TRUE) {data.type <- sheet}
    TYPES <- c("community", "environmental", "stacked")
    data.type <- match.arg(data.type, TYPES)
    if (is.null(sheet) == TRUE) {sheet <- data.type}
    if (data.type == "stacked") {
        stackeddata <- sqlFetch(dataplace,sheet)  
        result <- makecommunitydataset(stackeddata,row=sitenames,column=column,value=value,factor=factor,level=level)
    }else{
         result <- sqlFetch(dataplace,sheet,rownames=sitenames)
    }
    close(dataplace)
    rownames(result) <- make.names(rownames(result),unique=T)
    if (cepnames == TRUE) {
        colnames(result) <- make.cepnames(colnames(result))
    }else{
        colnames(result) <- make.names(colnames(result),unique=T)
    }
    return(result)
}

}