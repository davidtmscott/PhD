#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  January 2020
# Description:runs mFPT_new on landscapes 
#   from source to target (s.t) and from target to source (t.s). 
#   also calculates variance


# conditional first passage time

#file_name = "s10.t10"

source("mFPT.functions.R")

for (set in list(#list(n.target = 1, n.source = 1),
  #list(n.target = 100, n.source = 1),
  #list(n.target = 10, n.source = 10),
  #list(n.target = 1, n.source = 100),
  #list(n.target = 100, n.source = 100),
  #list(n.target = 1, n.source = 1000)
  #list(n.target = 1, n.source = 10),
  list(n.target = 10, n.source = 1))){
  #list(n.target = 10, n.source = 100),
  #list(n.target = 100, n.source = 10)
  
  print(set$n.source)
  print(set$n.target)
  
  n.source = set$n.source
  n.target = set$n.target
  #n.source = 100
  #n.target = 100
  file_name <- paste0("S", n.source, ".T", n.target)
  
  ##

  for (landscape.type in c("uniform", "clumpy", "corridor.1", "corridor.2")){
  #for (landscape.type in c("corridor.2")){
    print(landscape.type)
    
    rm(all.landscapes)
    
    load(file=paste0("../Data/", landscape.type,".landscapes.RData"))
    # load results lists as all.results
    load(file=paste0("../Results/Models/",file_name,"/",file_name,".",
                     landscape.type,".models.RData"))
    
    #landscape.type="uniform"
    
    # generate W and calculate mFPT
    for (i in 1:length(all.landscapes)){
      a = all.landscapes[[i]]
      b = all.results[[i]]
      
      mFPT_new.s.t.m = c()
      mFPT_new.s.t.m2 = c()
      #mFPT_new.s.t.v = c()
      
      mFPT_new.t.s.m = c()
      mFPT_new.t.s.m2 = c()
      #mFPT_new.t.s.v = c()
      
      for (j in 1:length(a$landscapes)){
        # add sources and targets to landscape
        a$landscapes[[j]] <- source.target(a$landscapes[[j]], 
                                           b$coords$source_coords, 
                                           b$coords$target_coords)
        W <- with(b$pars, trans_probs(a$landscapes[[j]], dispersal, R, total.patches, extent))
        sources <- b$coords$sources
        targets <- b$coords$targets
        internals <- as.numeric(colnames(W[,-c(sources, targets)]))
        
        mFPT_new_st <- new_mFPT(W, sources, internals, targets)
        mFPT_new.s.t.m <- c(mFPT_new.s.t.m, mean(mFPT_new_st$T1))
        mFPT_new.s.t.m2 <- c(mFPT_new.s.t.m2, mean(mFPT_new_st$T2))
        #mFPT_new.s.t.v <- c(mFPT_new.s.t.v, mean(mFPT_new_st$V))
        
        # swap target and source    
        mFPT_new_ts <- new_mFPT(W, targets, internals, sources)
        mFPT_new.t.s.m <- c(mFPT_new.t.s.m, mean(mFPT_new_ts$T1)) 
        mFPT_new.t.s.m2 <- c(mFPT_new.t.s.m2, mean(mFPT_new_ts$T2)) 
        #mFPT_new.t.s.v <- c(mFPT_new.t.s.v, mean(mFPT_new_ts$V))     
        
      }
      
      #mFPT_new.min.all <- pmin(mFPT_new.s.t.all, mFPT_new.t.s.all)
      #x <- as.numeric(names(mFPT_new.min.all))
      #mFPT_new.which.all <- replace(x, x>1, 2) 
      
      all.results[[i]]$mFPT_new.s.t <- mFPT_new.s.t.m
      all.results[[i]]$mFPT_new.s.t.m2 <- mFPT_new.s.t.m2
      #all.results[[i]]$mFPT_new.s.t.v <- mFPT_new.s.t.v
      
      all.results[[i]]$mFPT_new.t.s <- mFPT_new.t.s.m
      all.results[[i]]$mFPT_new.t.s.m2 <- mFPT_new.t.s.m2
      #all.results[[i]]$mFPT_new.t.s.v <- mFPT_new.t.s.v
      
      #all.results[[i]]$mFPT_new.min <- mFPT_new.min.all
      #all.results[[i]]$mFPT_new.which <- mFPT_new.which.all
      
    }
    save(all.results, file=paste0("../Results/Models/",file_name,"/",file_name,".",
                                  landscape.type,".models.RData"))
    print("Done")
    
  }
}
