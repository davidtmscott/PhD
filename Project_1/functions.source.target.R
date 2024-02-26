#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date: November 2022
# Description: function to bind source adn target to landscape


source.target <- function(x, source_coords, target_coords){
  x = as.data.frame(x)
  x = rbind(source_coords, x)
  x = rbind(x, target_coords)
  return(x)
}

# seq(0,10,length.out = 250)



