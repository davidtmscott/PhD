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
rm(list = ls())

n.source = 250
n.target = 250

if (n.source == 1){
    source_coords <- c(5, 0)
  } else{
    source_coords <- data.frame(x=seq(from = 0, to = 10, length.out = n.source),
                                y=rep(0, n.source))
} 


if (n.target == 1){
  target_coords <- c(5,10)
} else{
  target_coords <- data.frame(x=seq(from = 0, to = 10, length.out = n.target),
                              y=rep(10, n.target))
} 




for (landscape.type in c("uniform", "clumpy", "corridor.1", "corridor.2")){
  print(landscape.type)
  
  #landscape.type = "uniform"
  
  file_name <- paste0("S", n.source, ".T", n.target)
  
  # load data 
  # object saved as all.landscapes
  load(file=paste0("../Data/", landscape.type,".landscapes.RData"))
  
  all.results <- list()
  for (i in 1:length(all.landscapes)){
    #i = 1
    a = all.landscapes[[i]]
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
  save(all.results, file=paste0("../Results/Models/",file_name,"/",file_name,".",
                                landscape.type,".models.RData"))
  print("done")
}

  
###





