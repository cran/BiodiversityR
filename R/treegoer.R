`treegoer.widen` <- function(
    treegoer,
    species=unique(treegoer$species)[1:100],
    filter.vars=c("bio05", "bio14", "climaticMoistureIndex")
)
{
  
  treegoer.subset <- treegoer$species %in% species
  treegoer <- treegoer[treegoer.subset, ]  
  
  focal.treegoer <- filter.vars[1]
  focal.ranges <- treegoer[treegoer$var == focal.treegoer, 
                           c("species", "n", "MIN", "Q05", "QRT1", "QRT3", "Q95", "MAX")]
  names(focal.ranges)[2:8] <- paste0(focal.treegoer, "_", names(focal.ranges)[2:8])
  ranges.lookup <- focal.ranges
  
  for (i in 2:length(filter.vars)) {
    
    focal.treegoer <- filter.vars[i]
    focal.ranges <- treegoer[treegoer$var == focal.treegoer, 
                             c("species", "n", "MIN", "Q05", "QRT1", "QRT3", "Q95", "MAX")]
    names(focal.ranges)[2:8] <- paste0(focal.treegoer, "_", names(focal.ranges)[2:8])
    
    # check - note that data was missing for some explanatory variables especially soil 
    cat(paste(focal.treegoer, "- ranges for all species:", all.equal(focal.ranges$species, ranges.lookup$species), "\n"))
    
    ranges.lookup <- dplyr::left_join(ranges.lookup,
                                      focal.ranges,
                                      by="species")
    
#   ranges.lookup <- cbind(ranges.lookup, focal.ranges[, c(2:8)]) # used in the Rpub, can not handle 
  }  
  
  return(ranges.lookup)
  
}

`treegoer.filter` <- function(
    site.data, treegoer.wide, 
    filter.vars=c("bio05", "bio14", "climaticMoistureIndex"), 
    limit.vars=c("Q05", "Q95")) 
{
  filtered.data <- treegoer.wide
  for (f in 1:length(filter.vars)) {
    focal.var <- filter.vars[f]
    LL <- paste0(focal.var, "_", limit.vars[1])
    filtered.data <- filtered.data[filtered.data[, LL] <= as.numeric(site.data[, focal.var]), ]
    UL <- paste0(focal.var, "_", limit.vars[2])
    filtered.data <- filtered.data[filtered.data[, UL] >= as.numeric(site.data[, focal.var]), ]
  }
  return(filtered.data) # modify the function to return the list of the suitable species
  #  return(nrow(filtered.data)) # return number of species as in TreeGOER manuscript
}   

`treegoer.score` <- function(
    site.data, site.species=treegoer.wide$species,
    treegoer.wide, 
    filter.vars=c("bio05", "bio14", "climaticMoistureIndex")
) 
{
  
  site.species <- data.frame(species=site.species)
  filteri <- data.frame(treegoer.wide)
  
  site.species$climate.score <- rep(-1, nrow(site.species))
  site.species[site.species$species %in% filteri$species, "climate.score"] <- 0
  
  first.n <- paste0(filter.vars[1], "_n")
  
  site.species <- dplyr::left_join(site.species,
                                   filteri[ , c("species", first.n)],
                                   by="species")
  
  names(site.species)[which(names(site.species)==first.n)] <- "n.TreeGOER"
  
  filteri <- treegoer.filter(site.data=site.data,
                             treegoer.wide=filteri,
                             filter.vars=filter.vars,
                             limit.vars=c("MIN", "MAX"))
  
  site.species[site.species$species %in% filteri$species, "climate.score"] <- 1
  
  filteri <- treegoer.filter(site.data=site.data,
                             treegoer.wide=filteri,
                             filter.vars=filter.vars,
                             limit.vars=c("Q05", "Q95"))
  
  site.species[site.species$species %in% filteri$species, "climate.score"] <- 2
  
  filteri <- treegoer.filter(site.data=site.data,
                             treegoer.wide=filteri,
                             filter.vars=filter.vars,
                             limit.vars=c("QRT1", "QRT3"))
  
  site.species[site.species$species %in% filteri$species, "climate.score"] <- 3
  
  return(site.species)
} 
