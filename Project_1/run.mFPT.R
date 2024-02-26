#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  January 2020
# Description:

# file_name <- "S1.T1"
# landscape.type <- "uniform"

source("functions.mFPT.R")
source("functions.source.target.R")
source("functions.transition.rate.R")

for (set in list(
  list(n.source = 1, n.target = 1),
  list(n.source = 100, n.target = 1),
  list(n.source = 100, n.target = 100),
  list(n.source = 1000, n.target = 1),
  list(n.source = 1000, n.target = 1000))){
  
  print(set$n.source)
  print(set$n.target)
  #set = list(n.target = 1000, n.source = 1000)
  
  n.source = set$n.source
  n.target = set$n.target
  #n.source = 100
  #n.target = 100
  file_name <- paste0("S", n.source, ".T", n.target)
  
  ##
  
  # generate W and calculate condatis
  for (landscape.type in c("uniform", "clumpy", "corridor.1", "corridor.2")){ # 
  #landscape.type = "clumpy"
  #for (landscape.type in c("corridor.2")){ # 
    print(landscape.type)
    
    rm(all.landscapes)
    
    load(file=paste0("../Data/", landscape.type,".landscapes.RData"))
    # load results lists as all.results
    #load(file=paste0("../Results/Models/",file_name,"/",file_name,".",landscape.type,".models.RData"))
    load(file=paste0("../Results/Models/files_to_use/",file_name,"/",file_name,".",landscape.type,".models.RData"))
    
    # generate W and calculate mFPT
    for (i in 1:length(all.landscapes)){
      #i=34
      a = all.landscapes[[i]]
      b = all.results[[i]]
      
      mFPT.s.t = c()
      mFPT.s.t.m2 = c()
      
      mFPT.t.s = c()
      mFPT.t.s.m2 = c()
      
      for (j in 1:length(a$landscapes)){
        # add sources and targets to landscape
        #j=2
        a$landscapes[[j]] <- source.target(a$landscapes[[j]], 
                                           b$coords$source_coords, 
                                           b$coords$target_coords)
        W <- with(b$pars, trans_probs(a$landscapes[[j]], dispersal, R, total.patches, extent))
        sources <- b$coords$sources
        targets <- b$coords$targets
        
        # ensures cannot jump directly from the source to the target
        W[sources,targets] = 0
        
        mFPT <- mFPT_cont(W, targets) 
        mFPT.s.t <- c(mFPT.s.t, mean(mFPT[sources]))
        #mFPT.s.t.m2 <- c(mFPT.s.t.m2, mean(mFPT_cont.2(W, targets, mFPT)[sources]))
        
        # swap target and source    
        mFPT. <- mFPT_cont(W, sources)
        mFPT.t.s <- c(mFPT.t.s, mean(tail(mFPT., length(targets))))#sum(mFPT[targets] / length(targets))) 
        #mFPT2 <- mFPT_cont.2(W, sources, mFPT.)
        #mFPT.t.s.m2 <- c(mFPT.t.s.m2, mean(tail(mFPT2, length(targets))))#sum( mFPT_cont.2(W, sources, mFPT)[targets] / length(targets)))
      }
      all.results[[i]]$mFPT.s.t <- mFPT.s.t
      #all.results[[i]]$mFPT.s.t.m2 <- mFPT.s.t.m2
      
      all.results[[i]]$mFPT.t.s <- mFPT.t.s
      #all.results[[i]]$mFPT.t.s.m2 <- mFPT.t.s.m2
  
    }
    #save(all.results, file=paste0("../Results/Models/",file_name,"/",file_name,".",landscape.type,".models.RData"))
    save(all.results, file=paste0("../Results/Models/files_to_use/",file_name,"/",file_name,".",landscape.type,".models.RData"))
    print("Done")
  }  
 
}   
# # generate W and calculate mFPT
# for (i in 1:length(all.landscapes)){
#   a = all.landscapes[[i]]
#   mFPT = c()
#   for (j in 1:length(a$landscapes)){
#     W <- with(a$pars, trans_probs(a$landscapes[[j]], dispersal, R, n.patches, extent))
#     # can put try except methods here and record if error
#     mFPT <- c(mFPT, with(a$pars, sum(mFPT_cont(W, target)[source] / length(source))))
#   }
#   all.landscapes[[i]]$mFPT <- mFPT
# }