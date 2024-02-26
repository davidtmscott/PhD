#!/usr/bin/env python

"""

Code to run a single occurance of the simulation

used for testing


~~~~~~ environmental parameters

extent : size of landscape
n_years : duration of simulation 

temporal_G : rate of mean temp change over time
spatial_G : latitudinal temp gradient. 0 if no gradient

temporal_var & spatial_var : flags for variation in temp over space and time

# covariance functions (kernels) and n samples to generate in 1&2D gaussian processes
L : length scale for RBF
V : variance for RBF
gp_kernel_1D & gp_kernel_2D : kernels for fft

~~~~~~ biological parameters

density_dep : density dependence, controls fraction of population which surive after deaths each time step. 
                larger a, smaller fraction survive

lmbda : dispersal decay constant. larger lmbda, shorter dispersal distances
dispersal_frac : fraction of population that disperse each time step

population growth rate
T0 : temp when PGR = 0 and R = 1
GD : 
Rmax :



To do:
    choose GD and Rmax
    set up with __main__
    set up for parallel


In space the zero y index (latitude) represents the south and it increases latitudinally north. However, python displays 
arrays with the zero index at the north    
Changed lines in : 
    temp_spatial
    southern_immigration
    params_dict2df
    range_margin
    climate_margin
    initial_population
check results_dict2df

"""

# import packages
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

from funcs_Simulation import *

import pickle 
import cProfile
import pstats
import time

########################################################
# run simulation
########################################################

#df = pd.read_csv("../Data/S0.05_T0.02_sim_params.csv")
#df = pd.read_csv("../Data/inv_speed_0.5/inv_speed_0.5_G_sim_params.csv")


#2**10 
#2**11

extent_y = 2**13 #+ 2**11 
extent_x = 2**10
n_years = 10
wrap_size = 552 # can tolerate a lambda as low as 0.05

spatial_temp_range = 10
spatial_G = spatial_temp_range/(extent_y - (wrap_size*3)) # - 2**11

#spatial_G = 0

temporal_SD = [1] # 0.25, 0.5, 
SD = np.asarray(temporal_SD)
temporal_A = SD * np.sqrt(n_years)

# for no spatial or temporal variance, set V = 0
params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [n_years], # int(row["n_years"])  list(range(100, 10000, 100)
        
        "temporal_G" : [0.02], # 0.02
        "V_t" : [0],   "L_t" : [0], #54.7722557505166     # 1

        "spatial_G" : [spatial_G], #0.05
        "V_s" : [0],  "L_s" : [0], # 3238.17232401242   # 1

        "Rmax" : [10], #20
        "lambda" : [0.2012067021798425], #0.1001535502363233 #0.2012067021798425 # 0.13369555710396816 #  0.40940734064736645
        "dispersal_frac" : [0.75], #0.1
        
        "density_dep" : [0.1], 
        "T0" : [0], "GD" : [1], "S" : [0.01], # 0.05
        "tol" : [1e-10] } 

#print(params_values)
# create df of all param combos
params_df = params_dict2df(params_values)

#params_df.to_csv("../Data/n_years.csv", index=False)

# all param sets 
#margins, pop_dens, Rs, speeds, Ks = multi_param_sets(params_df)
params = params_df.iloc[[0]].to_dict("records")[0]


kernel, kernel_fft, params = make_kernel(params) # dispersal kernel"
southern_immigrants = southern_immigration_new(kernel, params) # from south of simulated landscape

del kernel 

data = {"kernel_fft" : kernel_fft, "southern_immigrants" : southern_immigrants, "wrap_size" : params["wrap_size"]}

t0 = time.time()
results = one_param_set(params, data)
t1 = time.time()


