#!/usr/bin/env python

"""
code to create param combinations to run simulations without spatial or temporal gradients

doing for jennys new experimental design

does not use the results of the analytical solution

"""


##########################

"""
calculate range boundary based on minimum Rmax value to be used

# for smallest Rmax
# calculate K0 just south of climate margin
# use this K to define the population range margin going forward
"""

import pandas as pd
import numpy as np
from funcs_Simulation import params_dict2df
from funcs_Temp import temp_spatial_lats, temp_spatial_grad
from funcs_PopDynamics import pop_growth_rate_new


extent_y = 2**13 #- 552
extent_x = 2**10

wrap_size = 552 # can tolerate a lambda as low as 0.05
spatial_temp_range = 10
spatial_G = spatial_temp_range/(extent_y - wrap_size *3)

params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [0], 
        
        "temporal_G" : [0], "V_t" : [0],   "L_t" : [10], 
        "spatial_G" : [spatial_G], "V_s" : [0],  "L_s" : [10], 

        "Rmax" : [1.1], # this is most important
        "lambda" : [0], "dispersal_frac" : [0], 
        
        "density_dep" : [0.1], "T0" : [0], "GD" : [1], "S" : [0.01], "tol" : [1e-10] } 

params_df = params_dict2df(params_values)


#params_df = pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results_NEW.csv')


Ks = []
for i in range(params_df.shape[0]):
        print(i)

        #params = params_df.to_dict(orient='records')[0]
        params = params_df.iloc[i].to_dict()

        # temperature lattice with initial temp
        lat_grid = temp_spatial_lats(params["extent_y"]) 
        temp_grid = temp_spatial_grad(params, lat_grid)
        # mean temp at southern border
        T2 = np.mean(temp_grid, 1)[0]
        B = 1 / (params["S"] * (T2 - params["T0"]))
        # get R and K for row just south of climate margin
        R = pop_growth_rate_new(temp_grid, params["T0"], B, params["Rmax"])[params["y0"]-1,1]
        K = (R - 1) / (params["density_dep"] * R) 
        Ks += [K]

(min(Ks))
# 0.014352285636780126
# I hard coded this into the simulation


##########################

"""

first get the invasion speed (at Rmax) for multiple param sets


"""

import pandas as pd
from funcs_Simulation import params_dict2df

extent_y = 2**13 #- 552
extent_x = 2**10
n_years = 2000
wrap_size = 552 # can tolerate a lambda as low as 0.05
spatial_temp_range = 10

#spatial_G = spatial_temp_range/(extent_y - wrap_size *3)
spatial_G = 0

# for no spatial or temporal variance, set V = 0
params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [n_years], 
        
        "temporal_G" : [0], "V_t" : [0],   "L_t" : [10], 

        "spatial_G" : [spatial_G], "V_s" : [0],  "L_s" : [10], 

        "Rmax" : [1.1, 1.5, 2.5, 5, 10, 20],
        "lambda" : [1.1571953034774167, 0.40940734064736645, 0.20120670217984254, 0.13369555710396816, 0.10015355023632337], 
        "dispersal_frac" : [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9], 
        
        "density_dep" : [0.1], "T0" : [0], "GD" : [1], "S" : [0.01], "tol" : [1e-10] } 

#print(params_values)
# create df of all param combos
params_df = params_dict2df(params_values)

params_df["status"] = "unstarted"
params_df["seed"] = list(range(params_df.shape[0]))

params_df.to_csv("../Data/Jenny_sims_Data/no_g_params.csv", index=False)


# pip install --force-reinstall -v "numpy==1.21.1" pip install --force-reinstall -v "numpy==1.21.1" pip install --force-reinstall -v "numpy==1.21.1" pip install --force-reinstall -v "numpy==1.21.1"


### process the output from sims 

"""

Calculate invasion speeds


"""

# ------------------------------
import numpy as np
import pickle
import pandas as pd 
import matplotlib.pyplot as plt
import os , os.path 

from funcs_Simulation import calc_params


DIR = '../Results/Jenny_sims_Results/no_g_params/'

names = [name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]

numbers = [int(name.replace("output","").replace(".pkl","")) for name in names]

counter = 0
#i = 78
for i in np.sort(numbers): 
        print(i)

        infilename = DIR + "output" + str(i) + ".pkl"

        with open(infilename, "rb") as output_file:
                output = pickle.load(output_file)
        params = output["params"]

        P_lim_D = np.asarray(output["limits"]["P_lim_D"])
        
        if i == 0:
                df = pd.DataFrame(params, index = [0])
                df["inv_speed"] = (P_lim_D[-1] - P_lim_D[-1000]) / 1000
        else:        
                df = df.append(params, ignore_index=True)
                df["inv_speed"].iloc[counter] = (P_lim_D[-1] - P_lim_D[-1000]) / 1000

        counter += 1


#df.to_csv("../Results/Jenny_sims_Results/inv_speeds.csv", index=False)



"""

recompute simulations with spatial and temporal gradients

"""

extent_y = 2**13 
extent_x = 2**10
n_years = 20000
wrap_size = 552 # can tolerate a lambda as low as 0.05

spatial_temp_range = 10
spatial_G = spatial_temp_range/(extent_y - wrap_size * 3)


#i = 0.01
df_all = df.iloc[:0,:].copy()
# select invasion speeds that keep with with each temporal gradient (climate speed)
for i in [0, 0.01, 0.02, 0.03, 0.04]:
        sub = df.loc[df["inv_speed"] >= (i / spatial_G)]
        sub=sub.assign(temporal_G = i) 
        df_all = pd.concat([df_all, sub])
        
df_all["spatial_G"] = spatial_G
df_all["clim_speed"] = df_all["temporal_G"] / df_all["spatial_G"] 

df_all = calc_params(df_all)

df_all["seed"] = range(len(df_all["spatial_G"]))
df_all["status"] = "unstarted"
df_all["n_years"] = n_years

df_all.to_csv("../Data/Jenny_sims_Data/S_T_g_params.csv", index=False)

# ~~~~~ edit 
#df_all = pd.read_csv("../Data/Jenny_sims_Data/S_T_g_params.csv")
#df_all2 = calc_params(df_all)
#np.unique(df_all["density_dep"])
#np.unique(df_all["N"])
#df_all.to_csv("../Data/Jenny_sims_Data/S_T_g_params_NEW.csv", index=False)
#########################################################################################


""""

process output and calculate lag and if it has settled or not


"""

# ------------------------------
import numpy as np
import pickle
import pandas as pd 
import matplotlib.pyplot as plt
import os , os.path 

from funcs_Simulation import calc_params


DIR = '../Results/Jenny_sims_Results/S_T_g_params_NEW_2/'

names = [name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]

numbers = [int(name.replace("output","").replace(".pkl","")) for name in names]

counter = 0
counter2 = 0
for i in np.sort(numbers): 
        print(i)
        #i=455

        infilename = DIR + "output" + str(i) + ".pkl"

        with open(infilename, "rb") as output_file:
                output = pickle.load(output_file)
        params = output["params"]

        P_lim_D = np.asarray(output["limits"]["P_lim_D"])
        C_lim = np.asarray(output["limits"]["C_lim_lin"])

        iterations = np.count_nonzero(P_lim_D)

        params["time_steps"] = iterations

        # if it settled within 10,000
        if iterations < 20000: 
                P_lim_D = P_lim_D[P_lim_D != 0]
                C_lim = C_lim[C_lim != 0]
                #len(P_lim_D)
                #len(C_lim)
                lag = P_lim_D[-1] - C_lim[-1]
        # if lag didnt settle
        else:
                # select those that had average rate of change < 0.001
                P_lim_D = P_lim_D[P_lim_D != 0]
                C_lim = C_lim[C_lim != 0]
                lag = P_lim_D - C_lim
                #plt.plot(range(len(lag)),lag);plt.show()
                counter2 += 1
                if (lag[-100] - lag[-1])/100 < 0.001:
                        lag = P_lim_D[-1] - C_lim[-1]
                        
                else:
                        lag = np.nan
        if i == 0:
                df = pd.DataFrame(params, index = [0])
                #df["G_inv_speed"] = (P_lim_D[-1] - P_lim_D[-100]) / 100
                df["lag"] = lag
        else:        
                df = df.append(params, ignore_index=True)
                #df["G_inv_speed"].iloc[counter] = (P_lim_D[-1] - P_lim_D[-100]) / 100
                df["lag"].iloc[counter] = lag

        counter += 1

counter2

# select those that did not settle 
df_s = df[df['lag'].isnull()]
df_s["n_years"] = 50000
df_s.to_csv('../Data/Jenny_sims_Data/S_T_g_params_2.csv',index=False)


df["group_id"] = df.groupby(["Rmax", "lambda"], as_index=False).grouper.group_info[0]  #, "dispersal_frac"

# remove that that did not settle
df = df[df['lag'].notna()]


df.to_csv('../Results/Jenny_sims_Results/S_T_g_results_NEW_2.csv',index=False)

###
df=pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results.csv')
df.loc[(df["Rmax"] == 20) & (df["dispersal_frac"] == 0.1) & (df["lambda"] == np.unique(df["lambda"])[1]),:]

df_old=pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results_old.csv')
df_old.loc[(df_old["Rmax"] == 20) & (df_old["dispersal_frac"] == 0.1) & (df_old["lambda"] == np.unique(df_old["lambda"])[1]),:]
##

"""
get all param combos which keep up with climate change

# get all unique groups of lambda and dispersal frac

get max Rmax for each group and create another x Rmax values from the max to zero

simulate

plot speed ~ Rmax, with lines for each lambda, dispersal_frac group

then try to interpolate



"""

dff = pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results.csv')

# remove temporal_G == 0
dff = dff[dff["temporal_G"] > 0]

dff["group_id"] = dff.groupby(["dispersal_frac", "lambda"], as_index=False).grouper.group_info[0]  #, "dispersal_frac"


a = np.linspace(20, 10, 15)[:-1] 
b = np.linspace(10, 5, 20)[:-1]
c = np.linspace(5,  2.5, 20)[:-1]
d = np.linspace(2.5, 1.5, 25)[:-1]
e = np.linspace(1.5, 1.1, 25)[:-1]
f = np.linspace(1.1, 1.01, 5)

np.linspace(20, 1, 200)

Rmaxs = np.concatenate((a, b, c, d, e, f))

for i in range(max(dff["group_id"])+1):
        #i=1
        sub = dff[dff["group_id"]==i].reset_index(drop=True)
        #sub[["Rmax","lambda", "dispersal_frac", "temporal_G"]]
        #print(min(sub["Rmax"]))
        sub = sub.iloc[[0]]
        sub = sub.loc[sub.index.repeat(len(Rmaxs))].reset_index(drop=True)
        sub["Rmax"] = Rmaxs
        if i == 0:
                df_df = sub.copy(deep=True)
        else:        
                df_df = df_df.append(sub, ignore_index=True)
                

df_df = df_df.drop(columns=['inv_speed', 'clim_speed', 'time_steps', 'lag'])

df_df["temporal_G"] = 0
df_df["spatial_G"] = 0
df_df["n_years"] = 2000
df_df["status"] = "unstarted"

df_df = calc_params(df_df)

df_df["N"]
df_df["Rmax"]

df_df.to_csv("../Data/Jenny_sims_Data/R_speeds.csv", index=False)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import numpy as np
import pickle
import pandas as pd 
import matplotlib.pyplot as plt
import os , os.path 

from funcs_Simulation import calc_params

from funcs_Temp import temp_spatial_grad, temp_spatial_lats
from funcs_PopDynamics import pop_growth_rate_new

DIR = '../Results/Jenny_sims_Results/R_speeds/'

names = [name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]

numbers = [int(name.replace("output","").replace(".pkl","")) for name in names]

counter = 0
#i = 78
for i in np.sort(numbers): 
        print(i)

        infilename = DIR + "output" + str(i) + ".pkl"

        with open(infilename, "rb") as output_file:
                output = pickle.load(output_file)
        params = output["params"]

        P_lim_D = np.asarray(output["limits"]["P_lim_D"])
        
        if i == 0:
                df = pd.DataFrame(params, index = [0])
                df["Rmax_inv_speed"] = (P_lim_D[-1] - P_lim_D[-1000]) / 1000
        else:        
                df = df.append(params, ignore_index=True)
                df["Rmax_inv_speed"].iloc[counter] = (P_lim_D[-1] - P_lim_D[-1000]) / 1000

        counter += 1


df.to_csv("../Results/Jenny_sims_Results/R_speeds_calc.csv", index= False)

"""

#df = pd.read_csv("../Results/Jenny_sims_Results/R_speeds_calc.csv")

#for i in np.unique(df["group_id"]):
i=1
print(i)
# select group 
s=df[df["group_id"] == i]
s.shape[0]
#plt.plot(s["Rmax"], s["Rmax_inv_speed"])
plt.scatter(s["Rmax"], s["Rmax_inv_speed"])

tck = interpolate.splrep(s["Rmax_inv_speed"], s["Rmax"])
plt.plot(tck[0], tck[1])
plt.show()

len(tck[0])
# interpolate R that matches clim speed 
#R_local = np.interp(row["clim_speed"], np.flip(s["Rmax_inv_speed"]), np.flip(s["Rmax"]))

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline



# Generate some noisy data
x = s["Rmax_inv_speed"]
y = s["Rmax"]
x_new = np.linspace(min(x), max(x), 1000)


from scipy.interpolate import CubicSpline

# Perform smoothed interpolation
spline = UnivariateSpline(x, y, s=0.01)  # 's' controls smoothness, smaller 's' means smoother curve
uni_y_smooth = spline(x_new)

spline(6)


# Perform cubic spline interpolation
cs = CubicSpline(x,y)
cb_y_smooth = cs(x_new)
cs(6)


tck = interpolate.splrep(s["Rmax_inv_speed"],s["Rmax"])
x_new = np.linspace(min(s["Rmax_inv_speed"]), max(s["Rmax_inv_speed"]), 1000)
spl_y_smooth = interpolate.splev(x_new, tck)
interpolate.splev(6, tck)

s["Rmax_inv_speed"][15:17]
tck = interpolate.splrep(s["Rmax_inv_speed"][14:18],s["Rmax"][14:18])
x_new = np.linspace(min(s["Rmax_inv_speed"][14:18]), max(s["Rmax_inv_speed"][14:18]), 1000)
spl_y_smooth = interpolate.splev(x_new, tck)
interpolate.splev(6, tck)


#interpolate.splev.root(6, tck)
#interpolate.sproot((tck[0], tck[1] - 6, tck[2]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

for i in np.unique(df["group_id"]):
        s=df[df["group_id"] == i]
        s.sort_values("Rmax_inv_speed", inplace=True)
        x = s["Rmax_inv_speed"]
        y = s["Rmax"]
        x_new = np.linspace(min(x), max(x), 1000)


        # Perform cubic spline interpolation
        cs = CubicSpline(x,y)
        cb_y_smooth = cs(x_new)
        tck = interpolate.splrep(x,y)
        spl_y_smooth = interpolate.splev(x_new, tck)


        plt.scatter(x, y, label='Noisy Data', s = 5)
        plt.plot(x_new, spl_y_smooth, label='spl Spline Interpolation', color='red', linewidth=1)
        #plt.plot(x_new, uni_y_smooth, label='uni Spline Interpolation', color='green', linewidth=1)
        plt.plot(x_new, cb_y_smooth, label='Cubic Spline Interpolation', color='black', linewidth=1)
        plt.ylabel("Rmax")
        plt.xlabel("Invasion Speed")
        plt.legend()
        plt.show()


R = pop_growth_rate_new(temp_grid, row["T0"], B, row["Rmax"])[:,1]
R_new = np.linspace(min(R), max(R), 1000)

#R_speed_location = np.interp(R_local, np.flip(R), np.flip(np.asarray(range(len(R)))))
tck = interpolate.splrep(np.flip(R), np.flip(np.asarray(range(len(R)))))
R_speed_location = interpolate.splev(R_new, tck)
interpolate.splev(9.28623438, tck)

plt.scatter(R, np.asarray(range(len(R))), label='Noisy Data', s=1)
plt.plot(R_new, R_speed_location, label='spl Spline Interpolation', color='red', linewidth=1)
plt.show()

"""

###################################

from scipy import interpolate

df = pd.read_csv("../Results/Jenny_sims_Results/R_speeds_calc.csv")

dff = pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results.csv')

dff["R_speed_location"] = np.nan

dff_sub = dff[dff["temporal_G"] == 0].reset_index(drop=True).copy()
dff_sub["group_id"] = np.nan

# remove temporal_G == 0
dff = dff[dff["temporal_G"] > 0].reset_index(drop=True)

dff["group_id"] = dff.groupby(["dispersal_frac", "lambda"], as_index=False).grouper.group_info[0]  #, "dispersal_frac"

dff.sort_values("group_id", inplace = True)

dff.reset_index(inplace=True, drop=True)


df["group_id"] = df.groupby(["dispersal_frac", "lambda"], as_index=False).grouper.group_info[0]  #, "dispersal_frac"
df.sort_values("group_id", inplace = True)
df.reset_index(inplace=True, drop=True)


# ~~ change this to use cubicspline interpolation
# ~~ check if there is a difference between the pop range limit and the sudden drop in population
#  see if this diff explains the diff between pop range limit and predicted pop range limit 
# ~~ check if intervals between the maximum Rmax and 1 were correct to use (non evenly distributed)
# ~~ redo whole process with pop range limit defined by lowest K value at climatic limit.
# ratio between K and Kmax

"""
np.unique(df['lambda']==row['lambda'])

dff.groupby(['lambda','dispersal_frac']).size().reset_index().rename(columns={0:'count'})

dff.groupby(['lambda','dispersal_frac'])

np.unique(df['lambda'])


for i in np.unique(df["group_id"]):
        #i=1
        print(i)
        s=df[df["group_id"] == i]
        print(s[["lambda", "dispersal_frac"]].iloc[0])


for i in np.unique(dff["group_id"]):
        #i=1
        print(i)
        s=dff[dff["group_id"] == i]
        print(s[["lambda", "dispersal_frac"]].iloc[0])


round(df['lambda'][0],8)==row['lambda']
np.unique(dff["lambda"])
"""

# 238
for i in range(dff.shape[0]):
        #i=1
        print(i)
        # select row 
        row = dff.loc[i,:]
        # select group 
        #df[df[['lambda','dispersal_frac']].isin(row[['lambda','dispersal_frac']])]
        s=df[(round(df['lambda'],8)==round(row['lambda'],8))&(df['dispersal_frac']==row['dispersal_frac'])].dropna()
        s.sort_values("Rmax_inv_speed", inplace=True)
        #s=df[df["group_id"] == row["group_id"]]
        if s.empty:
                dff.loc[i,"R_speed_location"] = np.nan
        else:
                # interpolate R that matches clim speed 
                #R_local = np.interp(row["clim_speed"], np.flip(s["Rmax_inv_speed"]), np.flip(s["Rmax"]))
                tck = interpolate.splrep(s["Rmax_inv_speed"], s["Rmax"])
                R_local = interpolate.splev(row["clim_speed"], tck)

                lat_grid = temp_spatial_lats(row["extent_y"]) # temperature lattice with initial temp
                temp_grid = temp_spatial_grad(dict(row), lat_grid)

                T2 = np.mean(temp_grid, 1)[0]
                B = 1 / (row["S"] * (T2 - row["T0"]))

                R = pop_growth_rate_new(temp_grid, row["T0"], B, row["Rmax"])[:,1]

                #R_speed_location = np.interp(R_local, np.flip(R), np.flip(np.asarray(range(len(R)))))
                tck = interpolate.splrep(np.flip(R), np.flip(np.asarray(range(len(R)))))
                R_speed_location = interpolate.splev(R_local, tck)

                #plt.plot(np.asarray(range(len(R))),R);plt.show()
                dff.loc[i,"R_speed_location"] = R_speed_location - row["y0"] #R_speed_pop_margin = np.interp(R_local, np.flip(R), np.flip(np.asarray(range(len(R)))))

dff_merged = pd.concat([dff_sub, dff], ignore_index=True, sort=False)


# interpolate all values or just a range around where I am interested?
# Rmax values chosen



# get actual clim lim
# check B

dff_merged.to_csv("../Results/Jenny_sims_Results/predicted_lags.csv", index= False)

"""
plt.scatter(dff["lag"], dff["R_speed_location"])
plt.xlim(-200,0)
plt.ylabel("Predicted lag")
plt.xlabel("Measured lag")
plt.show()




extent_y = 2**13 #- 552
wrap_size = 552 # can tolerate a lambda as low as 0.05
spatial_temp_range = 10
spatial_G = spatial_temp_range/(extent_y - wrap_size *3)

temporal_Gs = np.array([0.01, 0.02, 0.03, 0.04])

clim_speeds = temporal_Gs / spatial_G

for i in np.unique(df["group_id"]):
        #i=127
        s=df[df["group_id"] == i]
        print(i)
        #print(s.shape)

        plt.plot(s["Rmax"][:], s["Rmax_inv_speed"][:])
        plt.scatter(s["Rmax"][:], s["Rmax_inv_speed"][:])
        plt.show()


        s[["Rmax","Rmax_inv_speed"]][-50:-15]
        np.interp(17.99, np.flip(s["Rmax_inv_speed"]), np.flip(s["Rmax"]))

        #np.interp(17.99, np.flip(s["Rmax"]), np.flip(s["Rmax_inv_speed"]))


"""

# get climate speeds
# for each group get the R value which would give the same invasion speed as each climate speed
# for each param set get temp grid, then R grid. 
# find location where R value is










##########################################################################
##########################################################################
##########################################################################
##########################################################################


#######################################################
# data for plot of local invasion speed as a function of temperature or space

import numpy as np
import pandas as pd

from funcs_Simulation import params_dict2df

extent_y = 2**13 #- 552
extent_x = 2**10
n_years = 2000
wrap_size = 552 # can tolerate a lambda as low as 0.05


spatial_temp_range = 10
spatial_G = spatial_temp_range/(extent_y - wrap_size *3)

# for no spatial or temporal variance, set V = 0
params_values = {"extent_x" : [extent_x], "extent_y" : [extent_y], "n_years" : [n_years], # int(row["n_years"])  list(range(100, 10000, 100)
        
        "temporal_G" : [0], 
        "V_t" : [0],   "L_t" : [10], 

        "spatial_G" : [spatial_G], 
        "V_s" : [0],  "L_s" : [10], 

        "Rmax" : [20],
        "lambda" : [0.13369555710396816], 
        "dispersal_frac" : [0.1], 
        
        "density_dep" : [0.1], 
        "T0" : [0], "GD" : [1], "S" : [0.01], 
        "tol" : [1e-10] } 

#print(params_values)
# create df of all param combos
params_df = params_dict2df(params_values)

params_df["status"] = "unstarted"
params_df["seed"] = list(range(params_df.shape[0]))

from funcs_Simulation import create_matrices
from funcs_Dispersal import make_kernel
from funcs_SouthernDispersal import southern_immigration_new
import pyfftw
from funcs_PopDynamics import pop_growth_rate_new

params = params_df.iloc[[0]].to_dict("records")[0]
kernel, kernel_fft, params = make_kernel(params) # dispersal kernel"
southern_immigrants = southern_immigration_new(kernel, params) # from south of simulated landscape

del kernel 

data = {"kernel_fft" : kernel_fft, "southern_immigrants" : southern_immigrants, "wrap_size" : params["wrap_size"]}
b = pyfftw.empty_aligned((params["extent_y"], params["extent_x"] //2 +1), dtype='complex128', n=pyfftw.simd_alignment) 
params, data = create_matrices(params, data, b)

data.keys()
data["temp_grid"]

R = pop_growth_rate_new(data['temp_grid'], params["T0"], params["B"], params["Rmax"])

import matplotlib.pyplot as plt
plt.plot(np.array(range(len(R[:,0]))), R[:,0]);plt.show()


params_df = params_df.loc[params_df.index.repeat(8192)].reset_index(drop=True)

params_df["Rmax"] = R[:,0]
params_df["Rmax_original"] = 20
params_df["spatial_G"] = 0

params_df["temperature"] = data["temp_grid"][:,0]

params_df = params_df.loc[params["y0"]-1000:params["y0"],:]


params_df.to_csv("../Data/Jenny_sims_Data/no_g_params_PLOT_DATA.csv", index=False)


R[:,0][params["y0"] - 722]




# get speed at lag from results 

import numpy as np
import pickle
import pandas as pd 
import matplotlib.pyplot as plt
import os , os.path 

from funcs_Simulation import calc_params


DIR = '../Results/Jenny_sims_Results/no_g_params_PLOT_DATA/'

names = [name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]

numbers = [int(name.replace("output","").replace(".pkl","")) for name in names]

counter = 0
#i = 78
for i in np.sort(numbers): 
        print(i)

        infilename = DIR + "output" + str(i) + ".pkl"

        with open(infilename, "rb") as output_file:
                output = pickle.load(output_file)
        params = output["params"]

        P_lim_D = np.asarray(output["limits"]["P_lim_D"])
        
        if i == 0:
                df = pd.DataFrame(params, index = [0])
                df["inv_speed"] = (P_lim_D[-1] - P_lim_D[-1000]) / 1000
        else:        
                df = df.append(params, ignore_index=True)
                df["inv_speed"].iloc[counter] = (P_lim_D[-1] - P_lim_D[-1000]) / 1000

        counter += 1

import matplotlib.pyplot as plt

df["spatial_G"][0] 

df_all = pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results.csv')

sub = df_all.loc[(df_all['temporal_G'] > 0) & (df_all['Rmax'] == 20) & (df_all['dispersal_frac'] == 0.1)  & (df_all['lambda'] == np.unique(df_all["lambda"])[1])].reset_index()

np.flip(-np.asarray(range(len(df["inv_speed"]))))

plt.plot(-np.flip(range(len(df["inv_speed"]))),df["inv_speed"])

plt.vlines(sub["lag"][0],0, sub["clim_speed"][0])
plt.vlines(sub["lag"][1],0, sub["clim_speed"][1])
plt.vlines(sub["lag"][2],0, sub["clim_speed"][2])
plt.vlines(sub["lag"][3],0, sub["clim_speed"][3])
plt.hlines(sub["clim_speed"][0],sub["lag"][0],-1000)
plt.hlines(sub["clim_speed"][1],sub["lag"][1],-1000)
plt.hlines(sub["clim_speed"][2],sub["lag"][2],-1000)
plt.hlines(sub["clim_speed"][3],sub["lag"][3],-1000)
plt.xlabel("space")
plt.ylabel("Speed")
plt.show()

plt.plot(df["Rmax"][950:],df["inv_speed"][950:]);plt.show()







##########################################################################
##########################################################################
##########################################################################
##########################################################################




# ~~~~ getting the local invasion speed at the lag 

"""
calculate the speed of invasion at the lag
"""

import numpy as np 
import pandas as pd

from funcs_Simulation import create_matrices, params_dict2df, calc_params
from funcs_Dispersal import make_kernel
from funcs_SouthernDispersal import southern_immigration_new
import pyfftw
from funcs_PopDynamics import pop_growth_rate_new
from funcs_Temp import *

df = pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results.csv')

df2 = df.loc[df["temporal_G"] > 0,]

#df2.loc[df2["lag"] > 0,]

df2["R_lag"] = np.nan
df2["temp_lag"] = np.nan

# reset index
df2 = df2.reset_index(drop=True)

for i in range(df2.shape[0]):
        #i=0
        print(i)
        
        params = df2.loc[i,:]
        params = dict(params)

        lag_0 = int(np.floor(params["lag"]))
        lag_1 = lag_0 + 1

        lat_grid = temp_spatial_lats(params["extent_y"]) # temperature lattice with initial temp
        temp_grid = temp_spatial_grad(params, lat_grid)

        #data.keys()
        temp_0 = temp_grid[:,1][params["y0"]+lag_0]
        temp_1 = temp_grid[:,1][params["y0"]+lag_1]

        b = (lag_0 - lag_1) / (temp_0 - temp_1) 
        a = lag_1 - (b * temp_1)
        temp = (params["lag"] - a) / b

        R_lag = pop_growth_rate_new(temp, params["T0"], params["B"], params["Rmax"])

        df2.loc[i,"R_lag"] = R_lag
        df2.loc[i,"temp_lag"] = temp

        #df2.loc[i,:]

        # interpolate the actual temperature 
        # interpolate the actual R 
        # check if interpolated temperature gives interpolated R in pop_growth_rate_new
        # then get the speed of actual R
        # check if it is equal to climate speed




# change Rmax to new Rmax, store original 
df2["Rmax_orig"] = df2["Rmax"]
df2["Rmax"] = df2["R_lag"] 
# change spatial G to 0, store original
df2["spatial_G_orig"] = df2["spatial_G"]
df2["spatial_G"] = 0
# change temporal G to 0, store original
df2["temporal_G_orig"] = df2["temporal_G"]
df2["temporal_G"] = 0
# change n years
df2["n_years"] = 2000

df2 = calc_params(df2)

df2.to_csv('../Data/Jenny_sims_Data/lag_inv_speed.csv',index=False)

"""
calculate speed at lag
"""

# get speed at lag from results 

import numpy as np
import pickle
import pandas as pd 
import matplotlib.pyplot as plt
import os , os.path 

from funcs_Simulation import calc_params


DIR = '../Results/Jenny_sims_Results/lag_inv_speed/'

names = [name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]

numbers = [int(name.replace("output","").replace(".pkl","")) for name in names]

counter = 0
#i = 78
for i in np.sort(numbers): 
        print(i)

        infilename = DIR + "output" + str(i) + ".pkl"

        with open(infilename, "rb") as output_file:
                output = pickle.load(output_file)
        params = output["params"]

        P_lim_D = np.asarray(output["limits"]["P_lim_D"])
        
        if i == 0:
                df = pd.DataFrame(params, index = [0])
                df["lag_inv_speed"] = (P_lim_D[-1] - P_lim_D[-1000]) / 1000
        else:        
                df = df.append(params, ignore_index=True)
                df["lag_inv_speed"].iloc[counter] = (P_lim_D[-1] - P_lim_D[-1000]) / 1000

        counter += 1



# match invasion speed at lag back with its original data 
# positive lags wont have a speed

"""
combine lag speed with rest of data 
"""

df_all = pd.read_csv('../Results/Jenny_sims_Results/S_T_g_results.csv')


df.columns

df_sub = df[["seed", "R_lag", "temp_lag", "lag_inv_speed"]]

df2 = pd.merge(df_all, df_sub, on="seed", how="left")

df_all.loc[df_all["seed"]==466]["B"]
df.loc[df["seed"]==466]["B"]
df2.loc[df2["seed"]==466]["B"]


df2.to_csv('../Results/Jenny_sims_Results/S_T_g_results_lag_speed.csv', index = False)

# ~~~~~~~~~~~~~~~~~~~#
# plot them
for i in np.sort(numbers): 
        print(i)
        #i=267

        infilename = DIR + "output" + str(i) + ".pkl"

        with open(infilename, "rb") as output_file:
                output = pickle.load(output_file)
        params = output["params"]

        P_lim_D = np.asarray(output["limits"]["P_lim_D"])
        C_lim = np.asarray(output["limits"]["C_lim_lin"])
        
        if np.any(P_lim_D):
                P_lim_D = P_lim_D[P_lim_D != 0]
                C_lim = C_lim[C_lim != 0]

        lag = P_lim_D - C_lim
        #if len(lag) > 9000:
        plt.plot(range(len(lag)),lag);plt.show()



plt.scatter(df["inv_speed"] / df["clim_speed"], df["lag"])
plt.ylabel("lag")
plt.xlabel("Invasion speed / climate speed")
plt.show()

plt.scatter(df["clim_speed"] , df["lag"])
plt.ylabel("lag")
plt.xlabel("Invasion speed ")
plt.show()


i=34
infilename = DIR + "output" + str(i) + ".pkl"
with open(infilename, "rb") as output_file:
        output = pickle.load(output_file)




