#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  January 2020
# Description: runs and saves condatis to lists for each landscape

source("functions.condatis.R")

for (set in list(#list(n.target = 1, n.source = 1),
  #list(n.target = 1, n.source = 250),
  #list(n.target = 1, n.source = 750),
  list(n.target = 1, n.source = 500))){

  print(set$n.source)
  print(set$n.target)
  
  n.source = set$n.source
  n.target = set$n.target
  #n.source = 100
  #n.target = 100
  file_name <- paste0("S", n.source, ".T", n.target)
  
  ##
  #landscape.type = "clumpy"
  for (landscape.type in c("uniform", "clumpy", "corridor.1", "corridor.2")){
    print(landscape.type)
    load(file=paste0("../Data/", landscape.type,".landscapes.RData"))
    # load results lists as all.results
    load(file=paste0("../Results/Models/",file_name,"/",file_name,".",
                     landscape.type,".models.RData"))
    
    # generate W and calculate condatis
    for (i in 1:length(all.landscapes)){
      a = all.landscapes[[i]]
      b = all.results[[i]]
      
      condatis = c()
      for (j in 1:length(a$landscapes)){
        # add sources and targets to landscape
        a$landscapes[[j]] <- source.target(a$landscapes[[j]], 
                                           b$coords$source_coords, 
                                           b$coords$target_coords)
        W <- with(b$pars, trans_probs(a$landscapes[[j]], dispersal, R, total.patches, extent))
        # can put try except methods here and record if error
        condatis <- c(condatis, with(b$coords, condatis_s(W, targets, sources)))
      }
      all.results[[i]]$condatis <- condatis 
    }
    save(all.results, file=paste0("../Results/Models/",file_name,"/",file_name,".",
                                  landscape.type,".models.RData"))
  }
}
  