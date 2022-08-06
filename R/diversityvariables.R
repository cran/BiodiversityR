`diversityvariables` <-
function(x, y, digits=8){
    y$richness <- diversityresult(x, index='richness', method='each site', digits=digits)[,1]
    y$Shannon <- diversityresult(x, index='Shannon', method='each site', digits=digits)[,1]
    y$Simpson <- diversityresult(x, index='Simpson', method='each site', digits=digits)[,1]
    y$inverseSimpson <- diversityresult(x, index='inverseSimpson', method='each site', digits=digits)[,1]
    y$simpson.unb <- diversityresult(x, index='simpson.unb', method='each site', digits=digits)[,1]
    y$simpson.unb.inverse <- diversityresult(x, index='simpson.unb.inverse', method='each site', digits=digits)[,1]
    y$Logalpha <- diversityresult(x, index='Logalpha', method='each site', digits=digits)[,1]
    y$Berger <- diversityresult(x, index='Berger', method='each site', digits=digits)[,1]
    y$Jevenness <- diversityresult(x, index='Jevenness', method='each site', digits=digits)[,1]
    y$Eevenness <- diversityresult(x, index='Eevenness', method='each site', digits=digits)[,1]
    y$richness <- diversityresult(x, index='richness', method='each site', digits=digits)[,1]
    return(y)
}

