`import.from.Excel` <-
function(file=file.choose(),sheet="community",sitenames="sites",column="species",value="abundance",factor="",level="") {
    if (!require(RODBC)) {stop("Requires package RODBC")}
    dataplace <- odbcConnectExcel(file)
    SHEETS <- c("community", "environmental", "stacked")
    sheet <- match.arg(sheet, SHEETS)
    if (sheet=="stacked") {
        stackeddata <- sqlFetch(dataplace,sheet)  
        result <- makecommunitydataset(stackeddata,row=sitenames,column=column,value=value,factor=factor,level=level)
    }else{
         result <- sqlFetch(dataplace,sheet,rownames=sitenames)
    }
    close(dataplace)
    rownames(result) <- make.names(rownames(result),unique=T)
    colnames(result) <- make.names(colnames(result),unique=T)
    return(result)
}

