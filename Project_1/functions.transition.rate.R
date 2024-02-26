#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  December 2020
# Description: transition rate matrix

## ***** dispersal and lambda names should be swapped in accordance with the manuscript text

################################
# transition matrix 
################################

# takes positions and creates a transition rate matrix W
trans_probs <- function(landscape, dispersal, R, total.patches, extent){

  # get distances between each point 
  distances = dist(landscape, diag = T, upper = T)
  
  # set up negative exponential dispersal kernel
  #dispersal <- 2
  #cellside <- 1 
  #R <- 100
  
  lambda <- 2 / dispersal # was x / 2
  #norm <- R * lambda ^ 2 / 2 / pi * cellside ^ 4
  #k <- extent^2 * lambda^2 / (2 * pi * total.patches) # my old one
  
  # calculate transition rates based on distance calculated 
  #W <- norm * exp(-lambda * distances)
  #W <- R * k * exp(-lambda * distances) # my old one
  W <- R * (lambda^2) * exp(-lambda * distances)
  return(as.matrix(W))
}



# 
# extent = 10
# n.patches = 1000
# distances = seq(1:10)
# R = 40
# dispersal = c(0.5, 0.75, 1, 1.25)
# 
# for (i in dispersal){
#   print("stsart")
#   print(i)
#   
#   lambda <- 2 / i # was x / 2
#   
#   print(lambda)
#   #norm <- R * lambda ^ 2 / 2 / pi * cellside ^ 4
#   k <- extent^2 * lambda^2 / (2 * pi * n.patches)
#   
#   # calculate transition rates based on distance calculated 
#   #W <- norm * exp(-lambda * distances)
#   W <- R * k * exp(-lambda * distances)
#   
#   plot(distances, W, type = "l")
#   print("plotted")
# }



