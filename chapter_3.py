#!/usr/bin/env python

"""
code to create param combinations to run simulations WITH spatial or temporal variations

"""


"""

# biological parameters 

inv_speed/clim_speed = ~2

28.00859 /  13.072  = 2.142640

select Rmax, d, lambda

Rmax : 20
d : 0.1
lambda : 0.13369555710396816

temporal_G : 0.02

lag : -50.27299
time steps to settle (no var) : 184

lag no climate change : 9.43895

# temporal variation

select temporal amplitude and elength scale 

    - find the V that gives SD +- 0.5, 1, 1.5

- simulate without climate change (temporal_G) == 0

- simulate with climate change 

# spatial variation

select spatial amplitude and length scale 

- simulate without climate change 

- simulate with climate change 

# spatial AND temporal variation

- simulate without climate change 

- simulate with climate change

"""



import pandas as pd
import numpy as np
from funcs_Simulation import params_dict2df, calc_params

extent_y = 2**13 + 2**11 
extent_x = 2**10
n_years = 3000
wrap_size = 552 # can tolerate a lambda as low as 0.05
spatial_temp_range = 10

spatial_G = spatial_temp_range/(extent_y - (wrap_size *3) - 2**11)

#temporal_G = 0.02
temporal_G = [0, 0.01, 0.02, 0.03, 0.04]

Rmax =  20
lmbda = 0.13369555710396816
dispersal_frac = 0.1

# ~~~ code to check lags with smooth gradients
#df=pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results.csv')
#df.loc[(df["Rmax"] == 20) & (df["dispersal_frac"] == 0.1) & (df["lambda"] == np.unique(df["lambda"])[1]),:]


temporal_SD = [0.25, 0.5, 1]
temporal_L = [1, 3, 10, 30, 100]

spatial_SD = [0.25, 0.5, 1]
spatial_L = [1, 3, 10, 30, 100]


"""
Spatial variation WITHOUT climate change
"""


SD = np.asarray(spatial_SD)
spatial_A = SD * np.sqrt(extent_x * extent_y)


# for temporal variance
params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [n_years], # int(row["n_years"])  list(range(100, 10000, 100)
        
        "temporal_G" : temporal_G, 
        "V_t" : [0],   "L_t" : [0], 

        "spatial_G" : [spatial_G], 
        "V_s" : list(spatial_A),  "L_s" : spatial_L, 

        "Rmax" : [Rmax],
        "lambda" : [lmbda], 
        "dispersal_frac" : [dispersal_frac], 
        
        "density_dep" : [0.1], 
        "T0" : [0], "GD" : [1], "S" : [0.01], 
        "tol" : [1e-10] } 

#print(params_values)
# create df of all param combos
spatial_NO_cc = params_dict2df(params_values)

spatial_NO_cc = spatial_NO_cc.loc[np.repeat(spatial_NO_cc.index, 10)].reset_index(drop=True)

spatial_n = spatial_NO_cc.shape[0]
spatial_NO_cc["seed"] = range(spatial_n)#range(temporal_n, temporal_n +spatial_n)

#"""
#Spatial variation WITH climate change
#"""
#
#spatial_with_cc = spatial_NO_cc.copy()
#spatial_with_cc["temporal_G"] = temporal_G


"""
Temporal variation WITHOUT climate change
"""

SD = np.asarray(temporal_SD)
temporal_A = SD * np.sqrt(n_years)

# for temporal variance
params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [n_years], # int(row["n_years"])  list(range(100, 10000, 100)
        
        "temporal_G" : temporal_G, 
        "V_t" : list(temporal_A),   "L_t" : temporal_L, 

        "spatial_G" : [spatial_G], 
        "V_s" : [0],  "L_s" : [0], 

        "Rmax" : [Rmax],
        "lambda" : [lmbda], 
        "dispersal_frac" : [dispersal_frac], 
        
        "density_dep" : [0.1], 
        "T0" : [0], "GD" : [1], "S" : [0.01], 
        "tol" : [1e-10] } 

# create df of all param combos
temporal_NO_cc = params_dict2df(params_values)

temporal_NO_cc = temporal_NO_cc.loc[np.repeat(temporal_NO_cc.index, 10)].reset_index(drop=True)

temporal_n = temporal_NO_cc.shape[0]
temporal_NO_cc["seed"] = range(spatial_n, spatial_n +temporal_n)

#"""
#Temporal variation WITH climate change
#"""
#
#temporal_with_cc = temporal_NO_cc.copy()
#temporal_with_cc["temporal_G"] = temporal_G


"""
Spatial AND temporal variation WITHOUT climate change
"""


# for temporal variance
params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [n_years], # int(row["n_years"])  list(range(100, 10000, 100)
        
        "temporal_G" : [0], 
        "V_t" : list(temporal_A),   "L_t" : temporal_L, 

        "spatial_G" : [spatial_G], 
        "V_s" : list(spatial_A),  "L_s" : spatial_L, 

        "Rmax" : [Rmax],
        "lambda" : [lmbda], 
        "dispersal_frac" : [dispersal_frac], 
        
        "density_dep" : [0.1], 
        "T0" : [0], "GD" : [1], "S" : [0.01], 
        "tol" : [1e-10] } 

# create df of all param combos
temp_spatial_NO_cc = params_dict2df(params_values)

temp_spatial_NO_cc = temp_spatial_NO_cc.loc[np.repeat(temp_spatial_NO_cc.index, 10)].reset_index(drop=True)

temp_spatial_n = temp_spatial_NO_cc.shape[0]
temp_spatial_NO_cc["seed"] = range(temporal_n + spatial_n, temporal_n +spatial_n+temp_spatial_n)


"""
Spatial AND temporal variation WITH climate change
"""


temp_spatial_with_cc = temp_spatial_NO_cc.copy()
temp_spatial_with_cc["temporal_G"] = temporal_G

"""
combine them all
"""

#temporal_NO_cc
#temporal_with_cc

#spatial_NO_cc
#spatial_with_cc

#temp_spatial_NO_cc
#temp_spatial_with_cc

#df = pd.concat([temporal_with_cc, temporal_NO_cc, 
#                spatial_with_cc, spatial_NO_cc, 
#                temp_spatial_with_cc, temp_spatial_NO_cc])

df = pd.concat([spatial_NO_cc, temporal_NO_cc])

df["status"] = "unstarted"

#df = calc_params(df)

df["clim_speed"] = df["temporal_G"] / df["spatial_G"] 


df.reset_index(drop=True, inplace=True)

#df.to_csv("../Data/Jenny_sims_Data/S_T_var_params.csv", index=False)

######################################################################################

# process and plot

# new simulations 
#       create new data and upload 
#       create new directory 
#       update scripts 
#       run 


"""
Process output 
"""


import numpy as np
import pickle
import pandas as pd 
import matplotlib.pyplot as plt
import os , os.path 


DIR = '../Results/Jenny_sims_Results/S_T_var_params_NEW_2/'

names = [name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]

numbers = [int(name.replace("output","").replace(".pkl","")) for name in names]

numbers = np.sort(numbers)#[:360]

numbers = numbers[:750-42]
numbers = numbers[750-42:]

counter = 0

for i in np.sort(numbers): 
        print(i)
        #i=210

        infilename = DIR + "output" + str(i) + ".pkl"

        with open(infilename, "rb") as output_file:
                output = pickle.load(output_file)
        params = output["params"]

        P_lim_D = np.asarray(output["limits"]["P_lim_D"]) #0.01
        C_lim = np.asarray(output["limits"]["C_lim_lin"])

        #iterations = np.count_nonzero(P_lim_D)
        #params["time_steps"] = iterations

        P_lim_D = P_lim_D[P_lim_D != 0]
        C_lim = C_lim[C_lim != 0]

        lag = P_lim_D - C_lim
        
        if i == 0:
                df = pd.DataFrame(params, index = [0])
                #df["G_inv_speed"] = (P_lim_D[-1] - P_lim_D[-100]) / 100
                df["mean_lag"] = np.nanmean(lag[-2000:])
        else:        
                df = df.append(params, ignore_index=True)
                #df["G_inv_speed"].iloc[counter] = (P_lim_D[-1] - P_lim_D[-100]) / 100
                df["mean_lag"].iloc[counter] = np.nanmean(lag[-2000:])

        counter += 1


df.to_csv("../Results/Jenny_sims_Results/S_T_var_results_all.csv", index=False)


# calculate mean and SE acrosss replicates 
df_means = df.groupby(["temporal_G", "V_t", "L_t", "V_s", "L_s"], as_index=False)["mean_lag"].mean()
df_means["n_reps"] = df.groupby(["temporal_G", "V_t", "L_t", "V_s", "L_s"], as_index=False).count()["seed"]
df_means["lag_SE"] = df.groupby(["temporal_G", "V_t", "L_t", "V_s", "L_s"], as_index=False)["mean_lag"].agg(np.std, ddof=1)["mean_lag"] / np.sqrt(df_means["n_reps"])

df_means.to_csv("../Results/Jenny_sims_Results/S_T_var_results_GROUPED_all_NEW_2_temporal.csv", index=False)

##########

"""
Extra simulations 

"""

# focus on old pop range limit 
# want extra climate speeds 
1 * 10 * 3 * 3
# want extra autocorr lengths 
3 * 10 * 2 * 3



"""
Spatial variation WITH climate change (temporal gradient 0.01)
"""



spatial_SD = [0.25, 0.5, 1]
spatial_L = [1,10,100]

SD = np.asarray(spatial_SD)
spatial_A = SD * np.sqrt(extent_x * extent_y)


params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [n_years], # int(row["n_years"])  list(range(100, 10000, 100)
        
        "temporal_G" : [0.01], 
        "V_t" : [0],   "L_t" : [0], 

        "spatial_G" : [spatial_G], 
        "V_s" : list(spatial_A),  "L_s" : spatial_L, 

        "Rmax" : [Rmax],
        "lambda" : [lmbda], 
        "dispersal_frac" : [dispersal_frac], 
        
        "density_dep" : [0.1], 
        "T0" : [0], "GD" : [1], "S" : [0.01], 
        "tol" : [1e-10] } 

#print(params_values)
# create df of all param combos
spatial_cc_01 = params_dict2df(params_values)

spatial_cc_01 = spatial_cc_01.loc[np.repeat(spatial_cc_01.index, 10)].reset_index(drop=True)

#spatial_n = spatial_cc_01.shape[0]
#spatial_cc_01["seed"] = range(spatial_n)




"""
Spatial variation WITH climate change EXTRA SPATIAL AUTOLENGTHSCALES
"""



spatial_SD = [0.25, 0.5, 1]
spatial_L = [3, 30]

SD = np.asarray(spatial_SD)
spatial_A = SD * np.sqrt(extent_x * extent_y)


params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [n_years], # int(row["n_years"])  list(range(100, 10000, 100)
        
        "temporal_G" : [0, 0.01, 0.02], 
        "V_t" : [0],   "L_t" : [0], 

        "spatial_G" : [spatial_G], 
        "V_s" : list(spatial_A),  "L_s" : spatial_L, 

        "Rmax" : [Rmax],
        "lambda" : [lmbda], 
        "dispersal_frac" : [dispersal_frac], 
        
        "density_dep" : [0.1], 
        "T0" : [0], "GD" : [1], "S" : [0.01], 
        "tol" : [1e-10] } 

#print(params_values)
# create df of all param combos
spatial_cc_L = params_dict2df(params_values)

spatial_cc_L = spatial_cc_L.loc[np.repeat(spatial_cc_L.index, 10)].reset_index(drop=True)

#spatial_n = spatial_cc_L.shape[0]
#spatial_cc_L["seed"] = range(spatial_n)


dff = pd.concat([spatial_cc_01, spatial_cc_L])

dff["status"] = "unstarted"

#dff = calc_params(dff)

dff["seed"] = np.asarray(range(1980, 1980 +dff.shape[0]))

dff["clim_speed"] = dff["temporal_G"] / dff["spatial_G"] 

dff.reset_index(drop=True, inplace=True)

#dff.to_csv("../Data/Jenny_sims_Data/S_T_var_params_2.csv", index=False)




###############################################################################

"""
process output 
"""


import numpy as np
import pickle
import pandas as pd 
import matplotlib.pyplot as plt
import os , os.path 


DIR = '../Results/Jenny_sims_Results/S_T_var_params_2_old/'

names = [name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]

numbers = [int(name.replace("output","").replace(".pkl","")) for name in names]

numbers = np.sort(numbers)#[:360]

counter = 0

for i in np.sort(numbers): 
        print(i)
        #i=210

        infilename = DIR + "output" + str(i) + ".pkl"

        with open(infilename, "rb") as output_file:
                output = pickle.load(output_file)
        params = output["params"]

        P_lim_D = np.asarray(output["limits"]["P_lim_D"]) #0.01
        C_lim = np.asarray(output["limits"]["C_lim_lin"])

        #iterations = np.count_nonzero(P_lim_D)
        #params["time_steps"] = iterations

        P_lim_D = P_lim_D[P_lim_D != 0]
        C_lim = C_lim[C_lim != 0]

        lag = P_lim_D - C_lim
        
        if i == 0:
                df = pd.DataFrame(params, index = [0])
                #df["G_inv_speed"] = (P_lim_D[-1] - P_lim_D[-100]) / 100
                df["mean_lag"] = np.nanmean(lag[-2000:])
        else:        
                df = df.append(params, ignore_index=True)
                #df["G_inv_speed"].iloc[counter] = (P_lim_D[-1] - P_lim_D[-100]) / 100
                df["mean_lag"].iloc[counter] = np.nanmean(lag[-2000:])

        counter += 1


#df.to_csv("../Results/Jenny_sims_Results/S_T_var_results_all_2.csv", index=False)


# calculate mean and SE acrosss replicates 
df_means = df.groupby(["temporal_G", "V_t", "L_t", "V_s", "L_s"], as_index=False)["mean_lag"].mean()
df_means["n_reps"] = df.groupby(["temporal_G", "V_t", "L_t", "V_s", "L_s"], as_index=False).count()["seed"]
df_means["lag_SE"] = df.groupby(["temporal_G", "V_t", "L_t", "V_s", "L_s"], as_index=False)["mean_lag"].agg(np.std, ddof=1)["mean_lag"] / np.sqrt(df_means["n_reps"])

df_means.to_csv("../Results/Jenny_sims_Results/S_T_var_results_GROUPED_all_2_old.csv", index=False)

