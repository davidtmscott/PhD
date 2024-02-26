#!/usr/bin/env python

"""
analytical solution to invasion speed numerically computed for different values of A

in derivation, A is lambda
"""

# import packages
import numpy as np
import matplotlib.pyplot as plt
from funcs_Dispersal import kernel_neg_exp_new

# set param set 
# calc R at carrying capacity 
# calculate dispersal as a sum

#def disp(D, A, disp_lmbda,extent_y):
#    Y = np.asarray(range(-extent_y+1, extent_y))
#    K = disp_lmbda * np.exp(-disp_lmbda * abs(Y)) / 2
#    K[extent_y] = 0 # set centre to 0
#    immigrants = D**2 * np.sum(K * A**(-Y))
#    return immigrants


def disp(params, A, Ky):
    if Ky.all() == None:
        print("MAKING NEW KENEL")
        # kernel depends on x and y 
        Kxy = kernel_neg_exp_new(params["extent_x"], params["lambda"])
        # kernel depends on y only
        Ky = np.sum(Kxy, 1) 
    Y = np.asarray(range(-params["extent_y"]+1, params["extent_y"]))
    Ky_padded = np.pad(Ky, ((len(Y) - len(Ky)) // 2, (len(Y) - len(Ky)) // 2), mode='constant', constant_values=0)
    imm = np.exp(np.log(Ky_padded) - Y * np.log(A))  
    immigrants = np.sum(imm)
    return immigrants

def omega(params, A, Ky):
    D = params["dispersal_frac"]
    R = params["Rmax"]
    w = R - (R * D)  + (R * D * disp(params, A, Ky))
    return w 

def speed(A, params, Ky):
    #A=0.1
    print(A)
    c = np.log(omega(params, A, Ky)) / np.log(A)
    return abs(c)



##
from scipy import optimize


# set invasion speed
extent_y = 2**13#7640 #2**13 - 552
extent_x = 2**10
wrap_size = 552 # can tolerate a lambda as low as 0.05
spatial_temp_range = 10
spatial_G = spatial_temp_range/(extent_y - wrap_size *3)


params_values = {"extent_x" : extent_x, "extent_y" : extent_y, "n_years" : 2000, # an_params["extent_y"] 

            "temporal_G" : 0, 
            "V_t" : 0,   "L_t" : 1, 

            "spatial_G" : 0, 
            "V_s" : 0,  "L_s" : 1, 
            
            "Rmax" : 5,
            "lambda" : 0.1001535502363233, # 0.2012067021798425 # 0.1001535502363233
            "dispersal_frac" : 0.01,
            
            "density_dep" : 0.1, 
            "T0" : 0, "GD" : 1, "S" : 0.05, # 0.05
            "tol" : 1e-10} 

# calc rest of params
params = calc_params(params_values)

Kernel = kernel_neg_exp_new(params["extent_x"], params["lambda"])
Kernel = cut_kernel(Kernel)

an_speed = optimize.minimize_scalar(speed, args=(params, Kernel))["fun"]

# lambda : 0.2012067021798425
# sim: 8.209908519007508
# an: 8.898040871687448

# lambda : 0.1001535502363233
# sim: 16.351905244444367
# an : 17
