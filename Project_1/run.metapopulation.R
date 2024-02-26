#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  January 2020
# Description: Generate W and run metapopulation simulation on landscapes
#   save in list

source("functions.metapopulation.R")
#n.target = 1
#n.source = 750

SE <- function(x) sd(x) / sqrt(length(x))

for (set in list(#list(n.target = 1, n.source = 1),
  #list(n.target = 1, n.source = 250),
  #list(n.target = 1, n.source = 750),
  list(n.target = 1, n.source = 750))){
  
  print(set$n.source)
  print(set$n.target)
  
  n.source = set$n.source
  n.target = set$n.target
  #n.source = 100
  #n.target = 100
  file_name <- paste0("S", n.source, ".T", n.target)
  
  ##
  
  #for (landscape.type in c("uniform", "clumpy", "corridor.1", "corridor.2")){
  for (landscape.type in c("clumpy")){
    # landscape.type = "clumpy"
    
    print(landscape.type)
    
    # load data 
    # object saved as all.landscapes
    load(file=paste0("../Data/", landscape.type,".landscapes.RData"))
    # load results lists as all.results
    load(file=paste0("../Results/Models/",file_name,"/",file_name,".",
                     landscape.type,".models.RData"))
    #i=1
    # generate W and run metapopulation simulations 
    for (i in 1:length(all.landscapes)){
      a = all.landscapes[[i]]
      b = all.results[[i]]
        
      metapop.mean = c()
      metapop.SE = c()
      metapop.var = c()
      metapop.sd = c()
      #j=1
      for (j in 1:length(a$landscapes)){
        # add sources and targets to landscape
        a$landscapes[[j]] <- source.target(a$landscapes[[j]], 
                                           b$coords$source_coords, 
                                           b$coords$target_coords)
        # transition rate of landscape with sources and targets
        W <- with(b$pars, trans_probs(a$landscapes[[j]], dispersal, R, total.patches, extent))
        metapop.raw = c()
        occupancy.frac = c()
        occupancy.total = c()
        for (k in 1:b$pars$n.metapops){
          mp <- with(b$pars, metapop.sim(total.patches, W, b$coords$sources, 
                                        b$coords$targets, n.source, extinction.rate))
          metapop.raw <- c(metapop.raw, mp$time)
          occupancy.frac <- c(occupancy.frac, mp$occupancy_frac)
          occupancy.total <- c(occupancy.total, mp$occupancy)
        }
        all.results[[i]]$metapop.raw[[j]] <- metapop.raw
        all.results[[i]]$occupancy.frac[[j]] <- occupancy.frac
        all.results[[i]]$occupancy.total[[j]] <- occupancy.total
        # calculate mean and se metapop
        metapop.mean = c(metapop.mean, mean(metapop.raw))
        metapop.SE = c(metapop.SE,  SE(metapop.raw))
        metapop.var = c(metapop.var, var(metapop.raw))
        metapop.sd = c(metapop.sd,  sd(metapop.raw))
      }
      all.results[[i]]$metapop.mean <- metapop.mean
      all.results[[i]]$metapop.SE <- metapop.SE
      all.results[[i]]$metapop.var <- metapop.var
      all.results[[i]]$metapop.sd <- metapop.sd
    }
    
    # save data 
    save(all.results, file=paste0("../Results/Models/",file_name,"/",file_name,".",
                                  landscape.type,".models.RData"))
    print("done")
  }

}

