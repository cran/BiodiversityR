`import.from.Access` <-
function(file=file.choose(),table="community",sitenames="sites",column="species",value="abundance",factor="",level="") {
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
    colnames(result) <- make.names(colnames(result),unique=T)
    return(result)
}

