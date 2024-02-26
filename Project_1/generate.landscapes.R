#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  October 2020
# Description: generate and save landscapes
#   sets paramets sets for landscapes
#   stores as lists of lists

rm(list = ls())

source("functions.landscape.R")

### generate and store landscapes for each landscape type and parameter set
#landscape.type = "uniform"
for (landscape.type in c("uniform", "clumpy", "corridor.1", "corridor.2")){
  print(landscape.type)
  
  # parameter values 
  vars = list(extent = 10, R = 1, extinction.rate = 0, n.metapops = 20,
              n.landscapes = 100, n.patches = c(200, 500, 1000, 2000), 
              dispersal = c(0.25, 0.5, 1, 2.5)) #0.7
  
  # landscape specific parameter values
  if (landscape.type == "uniform"){
    params= list(type = landscape.type )
  } else if (landscape.type == "clumpy"){
    params = list(type = landscape.type , n.foci = c(5, 10, 20), 
                  clump.length = c(0.25, 0.5, 1, 2)) #clump.frac = c(0.1, 0.25),
  } else if (landscape.type == "corridor.1"){
    params = list(type = landscape.type , n.foci = c(6, 12, 24), 
                  l.noise = c(0.1, 0.5, 1))
  } else if (landscape.type == "corridor.2"){
    params = list(type = landscape.type, 
                  n.foci = c(2, 5, 10), # was 5, 10, 20,
                  l.smooth = c(0.1, 1, 10),
                  n.trials.scaler = c(1.5, 3, 6))
  }
  vars = c(vars, params)
  
  # generate combinations
  vars = expand.grid(vars)
  #vars$target = vars$n.patches # update target
  
  if (landscape.type == "corridor.2"){
    vars$n.trials = vars$n.patches * vars$n.trials.scaler # change this
  }
  
  # generate landscapes for each set of parameter values
  all.landscapes <- list()
  for (i in 1:nrow(vars)){
    #i=1
    landscape.pars = Map(unlist, split.default(vars[i,], names(vars[i,])))
    landscape.pars = lapply(landscape.pars, unname)
    
    these.landscapes <- list()
    for (j in 1:landscape.pars$n.landscapes){
      if (landscape.type == "uniform"){
        this.landscape <- with(landscape.pars, 
                               data.frame(x = runif(n.patches, 0, extent), 
                                          y = runif(n.patches, 0, extent)))
      } else if (landscape.type == "clumpy"){
        this.landscape <- with(landscape.pars, 
                               simulate.clumpy.landscape(n.patches,extent, 
                                                         clump.length, n.foci))
      } else if (landscape.type == "corridor.1"){
        this.landscape <- with(landscape.pars, 
                               simulate.corridorey.landscape.1(n.patches,extent, 
                                                               n.foci, l.noise))
      } else if (landscape.type == "corridor.2"){
        this.landscape <- with(landscape.pars, 
                               simulate.corridorey.landscape.2(n.patches, extent, 
                                                               n.foci, l.smooth, 
                                                               n.trials))
      }
      #this.landscape <- source.target(this.landscape, source_coords, target_coords)
      these.landscapes <- c(these.landscapes, list(this.landscape))
    }
    the.landscapes <- list(pars=landscape.pars, landscapes=these.landscapes)
    all.landscapes <- c(all.landscapes, list(the.landscapes))
    
  }

  # save data 
  save(all.landscapes, file=paste0("../Data/", landscape.type,".landscapes.RData"))
  print("done")
}


#landscape.type = "corridor.2"
#load(file=paste0("../Data/", landscape.type,".landscapes.RData"))

# uniform = 16 
# clumpy = 192
# corridor.1 = 144
# corridor.2 = 432 




