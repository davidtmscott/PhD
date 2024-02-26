#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  February 2021
# Description: metapopulation function and test simulation


# distances between each patch is calculated 
# used to calculate dispersal probabilities between each patch 

# calculate Ti
# set parameter for number of dispersers from each patch
# simulation starts with one source patch occupied
# set time = 0
# at each time step: (gillespie algorithm)
#   calculate rate of each patch changing its state 
#   empty patch : sum of colonzation rates from all other occupied patched 
#   occupied patches : its rate of extinction
#   sum rates to see how short timestep will be
#   calculate Pij : rate (206 or 207 depending if empty or not) / sum rates to see how short timestep will be
#   choose which patch changes states based on Pij
#   generate T an exponetially dist random number with mean Ti and add to time taken
# simulation continuous until target site it colonized 
#    or maybe x% of target sites are occupied?

repeat{


################################################################################
# metapopulation simulation 

# W is a matrix, nxn
# occupancy is a vector size n: 0 for empty, 1 for occupied

meta.pop_step <- function(occupancy, extinction.rate, sources, 
                          col.rate.unconditional, W){
  #col.rate.unconditional <- occupancy %*% W  # a vector, size n
  # zero rates for occupied patches
  col.rate.conditional <- (1 - occupancy) * col.rate.unconditional
  # get extinction rate vector occupied patches (zero otherwise)
  extinction.rate.conditional <- occupancy * extinction.rate
  # combine 2 vectors so each patch has a rate
  transition.rate <- col.rate.conditional + extinction.rate.conditional
  # set source transition rate to 0 - cannot change state or go extinct
  transition.rate[sources] <- 0
  total.rate <- sum(transition.rate)
  # convert to probs
  transition.prob <- transition.rate / total.rate
  
  # select patch to change state
  which.transition <- sample(1:length(occupancy), 1, prob = transition.prob)
  # change occupancy of patch 
  occupancy.new <- occupancy
  occupancy.new[which.transition] <- 1 - occupancy.new[which.transition]
  # record time taken
  time <- rexp(1, total.rate)
  # changed state 
  changed.state <- occupancy.new[which.transition] - occupancy[which.transition]
  col.rate.unconditional <- col.rate.unconditional + (changed.state * W[,which.transition])
  return(list(occupancy = occupancy.new, time = time, 
              col.rate.unconditional = col.rate.unconditional))
}


metapop.sim <- function(total.patches, W, sources, targets, n.source, extinction.rate){
  #occupancy <- c(1, rep(0, n.patches - 1))
  occupancy <- c(rep(1,n.source), rep(0, total.patches - n.source))
  if (n.source > 1){
    col.rate.unconditional = rowSums(W[,1:n.source]) 
  }else{
      col.rate.unconditional = W[,1] 
  } 
  time = 0
  #while (occupancy[targets] != 1){
  while (!(1 %in% occupancy[targets])){
    a <- meta.pop_step(occupancy, extinction.rate, sources, 
                       col.rate.unconditional, W)
    occupancy <- a$occupancy
    time <- time + a$time
    col.rate.unconditional <- a$col.rate.unconditional 
  }
  occupancy_frac = sum(occupancy) / total.patches
  return(list(occupancy = sum(occupancy), time = time, 
              occupancy_frac = occupancy_frac))
}



# add a tryCatch wrapper to metapop.sim
metapop.sim_TRY <- function(total.patches, W, sources, targets, n.source, extinction.rate) {
  out <- tryCatch(
    {
      metapop.sim(total.patches, W, sources, targets, n.source, extinction.rate)
    },
    error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(list(time = NA, occupancy = NA))
    },
    warning=function(cond) {
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}


#metapop.sim_TRY(total.patches, W, sources, targets, n.source, extinction.rate=0)


break

# dont want anything else below
stop()

# block testing code 

source("functions.transition.rate.R")

source_coords = c(0.0001, 0.0002) # for condatis only
source_coords2 = c(0.0005, 0.0005) # for condatis only
target_coords = c(0.9001, 0.9002) # for mFPT and condatis   

n.patches = 10
extent = 10
dispersal = 0.1
R = 40
total.patches = n.patches
n.source = 2
n.target = 1
extinction.rate = 0

# create landscape of x,y coordinate points
landscape = data.frame(x = runif(n.patches-3, 0, extent), 
                       y = runif(n.patches-3, 0, extent))
landscape = rbind(source_coords, landscape)
landscape = rbind(source_coords2, landscape)
landscape = rbind(landscape, target_coords)


W <- trans_probs(landscape, dispersal, R, total.patches, extent)

targets = 10
sources = c(1,2)


##

# make landscapes with two coordinates the same each time 
# set source and target sites
source_coords = c(0.0001, 0.0002) # for condatis only
target_coords = c(0.9001, 0.9002) # for mFPT and condatis                       

metapop.results <- data.frame(matrix(ncol = 11, nrow = 0))

namess <- c("extent", "n.patches", "lambda", "extinction.rate", "R", 
            "cellside", "time", "occupancy", "occupancy.frac",
            "target.reached", "wall.time")

colnames(metapop.results) <- namess

# nested loops
for (extent in c(1)){ #seq(1, 10, 1)){
  
  for(n.patches in c(1000)){ #seq(100, 1000, 10)){
    
    for (lambda in c(2)){ #seq(2, 100, 2)){
      
      for (extinction.rate in c(20)){ #seq(2, 200, 2)){
        
        for (R in seq(1, 40, 1)){
          print(R)
          for (cellside in c(1)){#seq(5, 50, 10)){
            for (i in 1:5){
              
              # create landscape of x,y coordinate points
              landscape = data.frame(x = runif(n.patches - 2, 0, extent), 
                                     y = runif(n.patches - 2, 0, extent))
              landscape = rbind(source_coords, landscape)
              landscape = rbind(landscape, target_coords)
              
              W <- trans_probs(landscape, lambda, R, cellside, n.patches, extent)
              
              # initialize occupancy vector 
              occupancy <- c(1, rep(0, n.patches - 1))
              col.rate.unconditional = W[,1] 
              
              source = 1
              target = n.patches 
              time = 0
              
              run_time <- system.time({
                while (time <= 100){
                  #while (occupancy[target] != 1){
                  
                  a <- meta.pop_step(occupancy, extinction.rate, source, 
                                     col.rate.unconditional)
                  occupancy <- a$occupancy
                  time <- time + a$time
                  col.rate.unconditional <- a$col.rate.unconditional 
                }
              })
              occupancy_frac <- sum(occupancy) / n.patches
              target.reached <- occupancy[target]
              metapop.results[nrow(metapop.results) + 1, ] <- list(extent, n.patches, 
                                                                   lambda, 
                                                                   extinction.rate, R, 
                                                                   cellside, time, 
                                                                   sum(occupancy), 
                                                                   occupancy_frac,
                                                                   target.reached,
                                                                   as.numeric(run_time[3]))
            }
          }
        }
      }
    }
  }
}


old.kernel <- read.csv("../Results/extent.10.R.1.40.ext.20.csv")

old.kernel.1.10 <- rbind(metapop.results, old.kernel)

write.csv(old.kernel.1.10, "../Results/norm.kernel.extent.1.10.csv", row.names = F)

compare.with.other.comp


metapop.results.mean <- aggregate(occupancy.frac ~ R + cellside, data=metapop.results, FUN=mean)

metapop.results.mean$col_ext_ratio = unique(metapop.results.mean$R / unique(metapop.results$extinction.rate))
plot(metapop.results.mean$col_ext_ratio, metapop.results.mean$occupancy.frac, 
     col = metapop.results.mean$cellside)

mean_occupancy_frac = with(metapop.results, tapply(occupancy_frac, R, mean))
col_ext_ratio = unique(metapop.results$R / metapop.results$extinction.rate)
plot(col_ext_ratio, mean_occupancy_frac)

}
