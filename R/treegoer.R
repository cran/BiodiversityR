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
                           c("species", "n", "MIN", "Q05", "QRT1", "MEDIAN", "QRT3", "Q95", "MAX")]
  names(focal.ranges)[2:9] <- paste0(focal.treegoer, "_", names(focal.ranges)[2:9])
  ranges.lookup <- focal.ranges
  
  for (i in 2:length(filter.vars)) {
    
    focal.treegoer <- filter.vars[i]
    focal.ranges <- treegoer[treegoer$var == focal.treegoer, 
                             c("species", "n", "MIN", "Q05", "QRT1", "MEDIAN", "QRT3", "Q95", "MAX")]
    names(focal.ranges)[2:9] <- paste0(focal.treegoer, "_", names(focal.ranges)[2:9])
    
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
    upper.only.vars = NULL,
    lower.only.vars = NULL,
    limit.vars=c("Q05", "Q95")) 
{
  filtered.data <- treegoer.wide
  for (f in 1:length(filter.vars)) {
    focal.var <- filter.vars[f]

    if ((focal.var %in% upper.only.vars) == FALSE) {
      LL <- paste0(focal.var, "_", limit.vars[1])
      filtered.data <- filtered.data[filtered.data[, LL] <= as.numeric(site.data[, focal.var]), ]
    }
  
    if ((focal.var %in% lower.only.vars) == FALSE) {
      UL <- paste0(focal.var, "_", limit.vars[2])
      filtered.data <- filtered.data[filtered.data[, UL] >= as.numeric(site.data[, focal.var]), ]
    }
  
  }
  
  return(filtered.data) # modify the function to return the list of the suitable species
  #  return(nrow(filtered.data)) # return number of species as in TreeGOER manuscript
}   

`treegoer.score` <- function(
    site.data, site.species=treegoer.wide$species,
    treegoer.wide, 
    filter.vars=c("bio05", "bio14", "climaticMoistureIndex"),
    upper.only.vars = NULL,
    lower.only.vars = NULL
) 
{
  
  site.species <- data.frame(species=site.species)
  filteri <- data.frame(treegoer.wide)
  
  site.species$climate.score <- rep(-1, nrow(site.species))
  
  # only species with data
  first.n <- paste0(filter.vars[1], "_n")
  filteri <- filteri[is.na(filteri[, first.n]) ==  FALSE, ]
  
  site.species[site.species$species %in% filteri$species, "climate.score"] <- 0
  
  site.species <- dplyr::left_join(site.species,
                                   filteri[ , c("species", first.n)],
                                   by="species")
  
  names(site.species)[which(names(site.species)==first.n)] <- "n.TreeGOER"

  filteri <- treegoer.filter(site.data=site.data,
                             treegoer.wide=filteri,
                             filter.vars=filter.vars,
                             upper.only.vars=upper.only.vars,
                             lower.only.vars=lower.only.vars,
                             limit.vars=c("MIN", "MAX"))
  
  site.species[site.species$species %in% filteri$species, "climate.score"] <- 0.5
    
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

`treegoer.position` <- function(
    site.data, 
    treegoer.wide, 
    focal.var="bio01"
) 
{
  focal.test <- paste0(focal.var, c("_MIN", "_MEDIAN", "_MAX"))
  focal.test <- c("species", focal.test)
  treegoer.wide <- data.frame(treegoer.wide)
  data1 <- treegoer.wide[, focal.test]
  names(data1) <- c("species", "MIN", "MEDIAN", "MAX")
  
  data1$position <-  as.numeric(site.data[, focal.var]) - data1$MEDIAN
  
  data1$position_score <- ifelse(
    data1$position >= 0,
    data1$position / (data1$MAX - data1$MEDIAN),
    data1$position / (data1$MEDIAN - data1$MIN)
  )
  
  data1$position_score <- round(data1$position_score, 2)
  out.names <- c("species", "position_score")
  data1 <- data1[, out.names]
  names(data1)[2] <- paste0(focal.var, "_position")
  
  return(data1)
}

`treegoer.map` <- function(
    map.rast, map.species=treegoer[1, "species"],
    treegoer, 
    filter.vars=c("bio05", "bio14", "climaticMoistureIndex"),
    upper.only.vars = NULL,
    lower.only.vars = NULL,
    verbose=FALSE
) 
{
  if(all(filter.vars %in% names(map.rast)) == FALSE) {stop("Not all filter.vars in the rast")}
  
  spec.raster <- as.numeric(map.rast[[1]] < Inf)
  spec.raster3 <- spec.raster
  spec.raster2 <- spec.raster
  spec.raster1 <- spec.raster
  spec.raster05 <- spec.raster
  
  
  for (i in 1:length(filter.vars)) {
    
    var.focal <- filter.vars[i]
    if (verbose == TRUE) {cat(paste0(var.focal, " - "))} 
    
    spec.targets <- data.frame(treegoer[treegoer$species==as.character(map.species) & treegoer$var==var.focal, ])
    
    # 1. score 3
    
    raster.mask <- as.numeric(map.rast[[var.focal]] < as.numeric(spec.targets["QRT1"]))
    spec.raster3 <- terra::mask(spec.raster3, mask=raster.mask, maskvalue=1, updatevalue=0)
    
    raster.mask <- as.numeric(map.rast[[var.focal]] > as.numeric(spec.targets["QRT3"]))
    spec.raster3 <- terra::mask(spec.raster3, mask=raster.mask, maskvalue=1, updatevalue=0)
    
    # 2. score 2
    
    raster.mask <- as.numeric(map.rast[[var.focal]] < as.numeric(spec.targets["Q05"]))
    spec.raster2 <- terra::mask(spec.raster2, mask=raster.mask, maskvalue=1, updatevalue=0)
    
    raster.mask <- as.numeric(map.rast[[var.focal]] > as.numeric(spec.targets["Q95"]))
    spec.raster2 <- terra::mask(spec.raster2, mask=raster.mask, maskvalue=1, updatevalue=0)
    
    # 3. score 1
    
    raster.mask <- as.numeric(map.rast[[var.focal]] < as.numeric(spec.targets["MIN"]))
    spec.raster1 <- terra::mask(spec.raster1, mask=raster.mask, maskvalue=1, updatevalue=0)
    
    raster.mask <- as.numeric(map.rast[[var.focal]] > as.numeric(spec.targets["MAX"]))
    spec.raster1 <- terra::mask(spec.raster1, mask=raster.mask, maskvalue=1, updatevalue=0)
    
    # 4. score 0.5
    
    if ((var.focal %in% upper.only.vars) == FALSE) {
      raster.mask <- as.numeric(map.rast[[var.focal]] < as.numeric(spec.targets["MIN"]))
      spec.raster05 <- terra::mask(spec.raster05, mask=raster.mask, maskvalue=1, updatevalue=0)
    }
    
    if ((var.focal %in% lower.only.vars) == FALSE) {
      raster.mask <- as.numeric(map.rast[[var.focal]] > as.numeric(spec.targets["MAX"]))
      spec.raster05 <- terra::mask(spec.raster05, mask=raster.mask, maskvalue=1, updatevalue=0)
    }
    
    
  } # vars
  
  spec.raster <- spec.raster3 + spec.raster2 + spec.raster1 + spec.raster05
  
  rclmat <- as.matrix(data.frame(orig=c(0, 1, 2, 3, 4),
                                 updat=c(0, 0.5, 1, 2, 3)))
  
  spec.raster2 <- terra::classify(spec.raster, rclmat)  
  names(spec.raster2) <- paste0(as.character(map.species), " scores")  
  
  
  return(spec.raster2) 
}  
  
