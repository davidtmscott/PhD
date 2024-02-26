#!/usr/bin/env python

"""
Function for setting up and computing simulations

E.g simulate_time_step is what calls each subprocess at each iteration of the simulation
"""

# import packages
import numpy as np
#import matplotlib.pyplot as plt
from itertools import product
import pandas as pd

# import functions
from funcs_Temp import *
from funcs_PopDynamics import * 
from funcs_Dispersal import *
from funcs_SouthernDispersal import * 
from funcs_CollectResults import *

import pyfftw




def calc_params(df):
    """ calculates N, y0 and intercept """
    # df = params_df.loc[0,:] 
    # calculate for each combination of params 
    #df["N"] = (df["Rmax"] - 1) / (df["density_dep"] * df["Rmax"]) # equilibrium density
    # -------------------------
    # for new pop range margin :
    df["N"] = 10
    df["density_dep"] = (df["Rmax"] - 1) / (df["N"] * df["Rmax"]) 
    # -------------------------
    # starting position of climate range margin
    #if df["extent_y"] == 8192:
    #    df["y0"] = df["extent_y"] - (552*3)
    #elif df["extent_y"] == 10240:
    #    df["y0"] = df["extent_y"] - (552*3) - 2**11
    df["y0"] = 6536
    # ensures initial temp at y0 == T0
    df["intercept"] = df["T0"] + (df["spatial_G"] * df["y0"]) 
    return df


def params_dict2df(params_values):
    """ get unique combination of all param values and make each a row of df """
    combinations = product(*(params_values[i] for i in params_values))
    df = pd.DataFrame(list(combinations), columns=list(params_values.keys()))
    try:
        df2 = calc_params(df)
        return df2
    except:
        return df


def create_matrices(params, data, b):
    """ create a dict of matrices and lists for temp, pop, dis kernel and s immigrants """
    params["wrap_size"] = data["wrap_size"]
    # temporal
    # ~~ added +1 to n_years
    temps = temp_temporal(0, params["n_years"]+1, params) # temp increase over time
    if params["V_t"] != 0:
        temps_var = temporal_variation_fft(params["n_years"]+1, params["V_t"], params["L_t"], RBF)
    else:
        temps_var = np.zeros(params["n_years"]+1)
    data['temps'] = temps + temps_var
    data['temps_change'] = np.diff(data['temps'])
    # spatial
    data['lat_grid'] = temp_spatial_lats(params["extent_y"]) # temperature lattice with initial temp
    data['temp_grid'] = temp_spatial_grad(params, data['lat_grid'])
   
    # calculate B
    if params["spatial_G"] == 0 and params["temporal_G"] == 0:
        T2 = data["temp_grid"]
    else:
        T2 = np.mean(data["temp_grid"], 1)[0]
    params["B"] = 1 / (params["S"] * (T2 - params["T0"]))

    if params["V_s"] != 0:
        data['temp_var_grid'] = spatial_variation_fft(params["extent_x"], params["extent_y"], params["V_s"], params["L_s"], RBF)
        data['temp_grid'] += data['temp_var_grid'] 
    else:
        data['temp_var_grid'] = np.zeros((params["extent_y"], params["extent_x"]))

    # calc mean and var of spatial and temporal var 
    #params["spatial_V_mean"] = np.mean(abs(data['temp_var_grid']))
    #params["spatial_V_var"] = np.var(data['temp_var_grid'])
    
    #params["temporal_V_mean"] = np.mean(abs(temps_var))
    #params["temporal_V_var"] = np.var(temps_var)
    
    # store as data 
    #data['kernel_fft'] = kernel_fft
    #data['southern_immigrants'] = southern_immigrants 

    # initial population
    #data["R"] = pop_growth_rate(data['temp_grid'], params["T0"], params["GD"], params["Rmax"])
    if params["spatial_G"] > 0:
        R = pop_growth_rate_new(data['temp_grid'], params["T0"], params["B"], params["Rmax"])   
    else:
        R = np.full_like(data['temp_var_grid'], params["Rmax"])
 
    data['pop_grid'] = pyfftw.empty_aligned((params["extent_y"], params["extent_x"]), dtype='float64', n=pyfftw.simd_alignment) 
    data["fft"] = pyfftw.FFTW(data['pop_grid'], b, axes=(0,1), direction="FFTW_FORWARD", flags=("FFTW_MEASURE", ), threads=1) #, flags=("FFTW_MEASURE", )
    # inverse 
    data["ifft"] = pyfftw.FFTW(b, data['pop_grid'], axes=(0,1), direction='FFTW_BACKWARD', flags=("FFTW_MEASURE", ),threads=1) #, flags=("FFTW_MEASURE", )
    
    data['pop_grid'][:] = initial_population(R, params) # lattice of initial population 
    return params, data


def simulate_time_step(params, data, j, results):
    """ order matters """
    #j = 0
    #T = data['temp_grid'] + data['temps'][j] # update temperature 
    data['temp_grid'] += data['temps_change'][j] # update temperature 

    if params["spatial_G"] > 0:
        #data["R"] = pop_growth_rate(T, params["T0"], params["GD"], params["Rmax"])
        R = pop_growth_rate_new(data['temp_grid'], params["T0"], params["B"], params["Rmax"])
    else:
        R = np.full_like(data['pop_grid'], params["Rmax"])

    #data['pop_grid'] = dispersal(data, params)
    dispersal(params, data)
    #data['pop_grid'] = births(data['pop_grid'], data["R"])
    #data['pop_grid'] = deaths(data['pop_grid'], params["density_dep"]) # change name to surivors
    births_deaths(data['pop_grid'], R, params["density_dep"])
    collect_results(results, params, data, R, j)
    #return data

#10    4.291    0.429    4.291    0.429 M:\PhD\project2\Code\funcs_PopDynamics.py:60(deaths)
# 150.881 seconds


def shift_landscape(params, data, n, j):
    # if it has moved x places shift landscape with it
    # check where pop limit is too not to leave it behind
    # adjust starting point of T0

    # update population 
    update_pop_grid(data["pop_grid"], n) #data['pop_grid'] = 
    # update temperature
    update_lat_grid(data["lat_grid"], n) #data['lat_grid'] = 
    update_temp_var(data["temp_var_grid"], n) #data['temp_var_grid'] =
    data['temp_grid'] = np.add(temp_spatial_grad(params, data['lat_grid']), data['temp_var_grid']) 
    # add temp
    data['temp_grid'] += data['temps'][j+1] 
    # ~~ can rewrite this so only change new rows, not rebuilding entire grid
    #return data


def change_in_lat(series, params, data):
    """ """
    # find location of T0 or bio range margin 
    T_lim = np.floor(series)
    # let it move x places first as a burn in
    # find lat of most southernly point (this changes)
    S_lim = data["lat_grid"][0]
    n = int(T_lim - S_lim - params["y0"])
    # if T0 has moved, shift landscape
    return n   


import time
#import os, psutil, memory_profiler



def simulate_model_expanding(params, data): 
    """ """
    
    results = create_results_dic(params)

    #process = psutil.Process()
    #print(process.memory_info().rss)  # in bytes 

    counter = []
    
    for j in range(params["n_years"]):
    #j = 0
    #while True:

        print(j)
        #j=1
        
        #st = time.time()

        simulate_time_step(params, data, j, results) #data = 
        #collect_results(results, params, data, j) #results = 

        # calc any change in lat
        try:
            n = change_in_lat(results["limits"]["C_lim_lin"][j], params, data)
        except:
            # if no temp gradient, use pop range limit to adjust landscape
            n = change_in_lat(results["limits"]["P_lim_D"][j], params, data)

        # if T0 moved NORTH of Y0 by x lats, shift landscape north n lat distance
        if n >= 1: 
            #print("shifted" + str(n))
            shift_landscape(params, data, n, j) #data = 
            n = 0 # reset counter

        if j > 100:
            P_lim_D = results["limits"]["P_lim_D"][results["limits"]["P_lim_D"] != 0][-100:]
            C_lim = results["limits"]["C_lim_lin"][results["limits"]["C_lim_lin"] != 0][-100:]
            lag = P_lim_D - C_lim
            #if all(abs(i) < 0.001 for i in np.diff(lag)): #0.0005
            #    break
            # average diff in lag
            if (lag[-100] - lag[-1])/100 < 0.001:
                break
            
        #j+=1
                    
        #process = psutil.Process()
        #results["limits"]["P_lim_k"][j] = process.memory_info().rss
        #print(process.memory_info().rss)  # in bytes 
        #et = time.time()
        #print(et - st)

        # save the number of iterations until it settles
    print(j)
    results["temps"].extend(data["temps"][:params["n_years"]]) # #len(r_dic["Nlat_pop_limit"])
    return results

  
#np.mean(data['pop_grid'],1)

#plt.plot(range(params["extent"]), np.mean(data['pop_grid'],1));plt.show()

#np.gradient(range(params["extent"]), np.mean(data['pop_grid'],1))

#np.mean(data['temp_grid'],1)

#np.mean(data['temp_grid'] + data['temps'][j],1)

# if n_lat t1 - n_lat t0 < x 
#results["limits"]["P_lim_D"][j] - results["limits"]["P_lim_D"][j-1]
# add 0 to list 
# if list has e.g. 100 zeros in a row, break loop

#(params["N"] - data["pop_grid"][0,0]) / params["N"]


def one_param_set(params, data):
    """ takes one param set and runs simulation, returning results"""
    b = pyfftw.empty_aligned((params["extent_y"], params["extent_x"] //2 +1), dtype='complex128', n=pyfftw.simd_alignment) 
    # create new data dict & store kernel and south immigration grid  
    params, data = create_matrices(params, data, b)
    # run simulation and store results
    results = simulate_model_expanding(params, data)
    results["params"] = params
    # turn dics into df for each param set 
    return results

def multi_param_sets(params_df):
    """ """
    all_margins = pd.DataFrame(); all_pop_dens = pd.DataFrame()
    all_Rs = pd.DataFrame(); all_Ks = pd.DataFrame(); speeds = np.array([])
    lmbdas = []
    #i = 0
    for i in range(params_df.shape[0]): # for each unique param set
        print(i)
        np.random.seed(3)
        params = params_df.iloc[[i]].to_dict("records")[0] # select param set
        # make new kernel and southern immigration grid only for new lambdas
        if params["lambda"] not in lmbdas:
            Kyx = make_kernel(params["extent_x"], params["extent_y"], params["lambda"]) # dispersal kernel
            southern_immigrants = southern_immigration(params) # from south of simulated landscape
            lmbdas += [params["lambda"]]

        # run sim on one param set and return results
        results = one_param_set(params, Kyx, southern_immigrants) 
        margins, pop_dens, Rs, speed, Ks = format_results(params, results)
        
        # store in main dfs for all param sets
        all_margins = pd.concat([all_margins, margins], axis=0, ignore_index=True)
        all_pop_dens = pd.concat([all_pop_dens, pop_dens], axis=0, ignore_index=True)
        all_Rs = pd.concat([all_Rs, Rs], axis=0, ignore_index=True)
        all_Ks = pd.concat([all_Ks, Ks], axis=0, ignore_index=True)
        speeds = np.append(speeds, speed)
    return all_margins, all_pop_dens, all_Rs, speeds, all_Ks



###

# old 


#j=0
#def simulate_model(params, r_dic, d_dic, R_dic):
#    """ """
#    data = create_matrices(params)
#    #for j in range(params["n_years"]):
#    for j in range(3):
#        # update population across landscape
#        data['pop_grid'], R = simulate_time_step(params, data, j)
#        r_dic, d_dic, R_dic = collect_results(r_dic, d_dic, R_dic, params, data, j, R)
#
#        #print("total_K - total_N")
#        #print(np.sum(carrying_capacity(R, params)) - np.sum(data["pop_grid"]))        
#
#    #r_dic["temps"].extend(data["temps"])
#    r_dic["temps"].extend(data["temps"][:len(r_dic["Nlat_pop_limit"])])
#    return r_dic, d_dic, R_dic 



