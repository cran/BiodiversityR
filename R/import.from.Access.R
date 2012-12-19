if(.Platform$OS.type == "windows") {

`import.from.Access` <-
function(file=file.choose(),table="community",sitenames="sites",column="species",value="abundance",factor="",level="",cepnames=FALSE) {
    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- odbcConnectAccess(file)
    TABLES <- c("community", "environmental", "stacked")
    table <- match.arg(table, TABLES)
    if (table=="stacked") {
        stackeddata <- sqlFetch(dataplace,table)  
        result <- makecommunitydataset(stackeddata,row=sitenames,column=column,value=value,factor=factor,level=level)
    }else{
         result <- sqlFetch(dataplace,table,rownames=sitenames)
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

`import.from.Access2007` <-
function(file=file.choose(),table="community",sitenames="sites",column="species",value="abundance",factor="",level="",cepnames=FALSE) {
    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- odbcConnectAccess2007(file)
    TABLES <- c("community", "environmental", "stacked")
    table <- match.arg(table, TABLES)
    if (table=="stacked") {
        stackeddata <- sqlFetch(dataplace,table)  
        result <- makecommunitydataset(stackeddata,row=sitenames,column=column,value=value,factor=factor,level=level)
    }else{
         result <- sqlFetch(dataplace,table,rownames=sitenames)
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