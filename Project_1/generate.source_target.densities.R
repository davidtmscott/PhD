#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  October 2022
# Description: 
#       create a new list structure, same as landscapes but without the actual 
#       landscapes. Add in source and target params and coords
#       use lists for rest of analysis to store model results. 
#       refer back to landscape lists for landscape coords


# set source and target sites

#file_name <- "Density"
file_name <- "Percentage.10"

for (landscape.type in c("uniform", "clumpy", "corridor.1", "corridor.2")){
  print(landscape.type)
  
  #landscape.type = "uniform"
  
  # load data 
  # object saved as all.landscapes
  load(file=paste0("../Data/", landscape.type,".landscapes.RData"))
  
  all.results <- list()
  for (i in 1:length(all.landscapes)){
    #i = 1
    a = all.landscapes[[i]]
    
    if (file_name == "Density"){
      
      area <- a$pars$extent^2
      patch.distances <- sqrt(area / a$pars$n.patches)
      
      coords = seq(0,10, by = patch.distances)
      n.source <- length(coords)
      n.target <- length(coords)
      
      source_coords <- data.frame(x=rep(0, n.source), y=coords)
      target_coords <- data.frame(x=rep(10, n.target), y=coords)
      
    }
    # percentage of n.patches
    else if(sub(".[1-9].*", "", file_name) == "Percentage"){
      
      # takes the percentage to use from file_name
      frac <- as.numeric(sub(".*Percentage.", "", file_name)) / 100
      n.source <- a$pars$n.patches * frac
      n.target <- a$pars$n.patches * frac
      coords = seq(0,10, length.out = n.source)
      
      source_coords <- data.frame(x=rep(0, n.source), y=coords)
      target_coords <- data.frame(x=rep(10, n.target), y=coords)
      
    } 
         
    # remove landscapes
    a$landscapes <- NULL 
    
    # update params 
    a$pars$n.target = n.target
    a$pars$n.source = n.source
    a$pars$total.patches = a$pars$n.patches + a$pars$n.source + a$pars$n.target
    
    # add source and target patch coordinates as new list elements
    a$coords$source_coords = source_coords
    a$coords$target_coords = target_coords
    
    # add source and target index positions
    a$coords$sources = seq(1, a$pars$n.source)
    a$coords$targets = seq(a$pars$total.patches-a$pars$n.target+1, a$pars$total.patches)
    
    this.element <- list(pars=a$pars, coords=a$coords)
    all.results <- c(all.results, list(this.element))
  }
  # save data as new file in Results
  save(all.results, file=paste0("../Results/",file_name,"/",file_name,".",
                                landscape.type,".models.RData"))
  print("done")
}
