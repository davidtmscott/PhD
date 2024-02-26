#!/usr/bin/env Rscript 

# Author: David Scott
# Contact: 
# Date:  October 2020
# Description: Condatis implementation

# Comparing mFPT and Condatis with a stochastic metapopulation model 
# using the same transition rates
# see 2011 hodgson paper


rm(list = ls()) # clears workspaces
graphics.off() # clears any images

################################################################################
# Parameters



################################################################################
# Functions 

source("functions.source.target.R")
source("functions.transition.rate.R")
source("functions.mFPT.R")
source("functions.condatis.R")
source("functions.metapopulation.R")
#source("functions.landscape.R")


################################################################################
# Generate landscapes

#source("generate.landscapes.R")


#source("generate.source_target.R")





################################################################################
# Compute metrics
  
source("run.metapopulation.R")
print("metapop done!!!")

#source("run.condatis.R")
#print("condatis done!!!!")

#source("run.mFPT.R")
#print("mFPT done!!!")

#source("run.mFPT_new.R")
#print("mFPT_new done!!!")
  

#Warning messages:
#  1: In grid.Call.graphics(C_upviewport, as.integer(n)) :
#  cannot pop the top-level viewport ('grid' and 'graphics' output mixed?)
#2: In UseMethod("depth") :
#  no applicable method for 'depth' applied to an object of class "NULL"


################################################################################
# Data cleaning
# create dataframes from lists in .RData files 

#source("clean.data.R")

################################################################################
# Compute GAMs and create LaTeX tables

#source("create.GAM.tables.R")


################################################################################
# Generate figures 

#source("create.figures.R")





