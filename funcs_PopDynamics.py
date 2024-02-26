#!/usr/bin/env python

###############################################################################################
# population dynamics
###############################################################################################

"""
Functions for pop dynamics 
Population growth rate which is a function of temperature
births, surival , carry capacity, initial density 
"""
# imports
import numpy as np
#import matplotlib.pyplot as plt

# functions 

def pop_growth_rate(t_grid, T0, GD, Rmax):
    """ Calculates pop growth rate (PGR) from temperature and converts to reproduction ratio R
    takes temperature grid 
    gamma (G) = how sensitive R is to temp, lower values give a less steep increase pop growth rate with temp, thus lower values of R
    delta (D) = duration of breeding, acts as a scaler for population growth rate
    beta (B) = when growth rate saturates at max, only used when T > T0. Higher = more saturation
    T0 = temp at which growth rate is 0
    Rmax = maximum growth rate, controls y value
    R = repro ratio. Biger R, more go south
    B limits number going south
    exp converts PGR to R """
    # calculate B from GD and Rmax 
    B = (GD) / np.log(Rmax) # GD steepness, saturates at exp(GD / B)
    print(B)
    PGR = t_grid - T0
    PGR_alt = PGR / (1 + B * PGR)
    PGR[t_grid > T0] = PGR_alt[t_grid > T0] # temps > T0
    R =  np.exp(GD * PGR) 
    return R


#def pop_growth_rate_new(t_grid, T0, B, Rmax):
#    GD = B * np.log(Rmax)
#    # calculate B from GD and Rmax 
#    #B = (GD) / np.log(Rmax) # GD steepness, saturates at exp(GD / B)
#    PGR = t_grid - T0
#    PGR_alt = PGR / (1 + B * PGR)
#    PGR[t_grid > T0] = PGR_alt[t_grid > T0] # temps > T0
#    R = np.exp(GD * PGR) # R
#    return R


# check if it is faster to pass the variable in as data["pop_grid"]
# do also for t_grid in p_g_r_new 
# and remove t_grid - T0

def pop_growth_rate_new(t_grid, T0, B, Rmax):
      t_above_T0 = t_grid > T0
      #PGR = np.where(t_above_T0, (t_grid - T0) / (1 + B * (t_grid - T0)), t_grid - T0)
      PGR = np.where(t_above_T0, (t_grid) / (1 + B * (t_grid)), t_grid)
      R = np.exp(B * np.log(Rmax) * PGR)
      return R

def births(p_grid, R):
    """ R - reproduction ratio """
    p_grid *= R 
    return p_grid

# 10    3.720    0.372    3.720    0.372 M:\PhD\project2\Code\funcs_PopDynamics.py:56(births)
def deaths(p_grid, a):
    """ calculate number of deaths using survival fraction 
    survival fraction based on density dependence (a) """
    #survival_frac = 1 / (1 + a * p_grid)
    #p_grid *= survival_frac 
    p_grid /= (1 + a * p_grid)
    return p_grid

#   10    2.331    0.233    2.331    0.233 M:\PhD\project2\Code\funcs_PopDynamics.py:62(deaths)

def births_deaths(p_grid, R, a):
    """ R - reproduction ratio """
    np.multiply(R, p_grid, out=p_grid) # births
    np.divide(p_grid, (1 + a * p_grid), out=p_grid) # deaths
    #return p_grid

def carrying_capacity(R, params):
    """ 
    check if this holds up with spatial variation
    """
    R_flat = R.ravel()
    # set a tol
    R_flat[np.where(R_flat < params["tol"])] = 0 # was 1e-20
    # array of zeros  
    K = np.zeros(R_flat.shape)
    # where R > 1
    R1_index = np.where(R_flat >= 1)[0] 
    R1 = R_flat[R1_index]
    K[R1_index] = (R1 - 1) / (params["density_dep"] * R1)
    #p_grid = np.reshape(p_grid, (-1,params["extent"])) # convert back to 2d
    K = np.reshape(K, (params["extent_y"], params["extent_x"])) # convert back to 2d
    # set negatives to zero
    #K[K < 0] = 0
    return K

def initial_population(R, params):
    """ """
    #t_grid=data["temp_grid"][0:20,:]
    #R = pop_growth_rate(t_grid, params["T0"], params["GD"], params["Rmax"])
    #R = pop_growth_rate_new(t_grid, params["T0"], params["B"], params["Rmax"])
    p_grid = carrying_capacity(R, params)
    #p_grid = (R - 1) / (params["density_dep"] * R)

    # if all values the same, i.e. no range margin
    if np.all(p_grid == np.ravel(p_grid)[0]):
        p_grid[int(params["y0"]):,:] = 0
        #p_grid = np.zeros((params["extent"],params["extent"]))
        #np.put(p_grid, 50, 1)
    return p_grid



#p_grid = data['pop_grid']
##########  plots 

def plot_R_PGR(x, y):
    """ """
    x = x.ravel()
    y = y.ravel()
    order = np.argsort(x)
    xs = np.array(x)[order]
    ys = np.array(y)[order]
    plt.plot(xs, ys); plt.show()

#plot_R_PGR(temp_grid, R)

def plot_R_params(T0, GD, Rmax):
    """for testing params"""
    #0.55, 6
    B = GD / np.log(Rmax)
    #B = 1
    #Rmax = np.exp(GD/B)
    #print(B)
    temps = np.linspace(-40,80)
    T_minus, T_plus = temps[temps<=T0], temps[temps>T0]
    R_minus = np.exp(GD * (T_minus - T0)) # < T0
    R_plus = np.exp((GD * (T_plus - T0)) / (1 + B * (T_plus - T0))) # > T0
    R = np.concatenate([R_minus, R_plus]) 
    #R = np.log(R)
    print(np.max(R))
    print(np.where(temps == 0))
    print(temps)
    plt.plot(temps,R); plt.show()

#plot_R_params(T0= 0, GD = 1, Rmax = 2)
#plot_R_params(T0 = 10, GD = 0.3, Rmax = 4)

#Rmax = 2
#GD = 1
#B = GD / np.log(Rmax)
#T = 0 # assuming T0 = 0
#N = (Rmax - 1) / (0.1 * Rmax)
#R = 0
#while (R < (Rmax- 0.0001)): 
##while (r < np.log(Rmax)): 
#    r = (GD * T) / (1 + B * T)
#    R = np.exp(r)
#    print(R)
#    #(R - 1) / (0.1 * R) 
#    T += 1
#
#np.log(Rmax)
######
#p1 = p_grid 

#R = pop_growth_rate(t_grid, params["T0"], params["GD"], params["Rmax"])
##p_grid = (R - 1) / (params["density_dep"] * R)
#
## unravel to 1D
#R_flat = R.ravel()
#R_flat[np.where(R_flat < 1e-20)] = 0
## array of zeros  
#p_grid = np.zeros(R_flat.shape)
## where R > 0
#R1_index = np.where(R_flat > 0)[0]
#R1 = R_flat[R1_index]
#len(R1)
#p_grid[R1_index] = (R1 - 1) / (params["density_dep"] * R1)
#p_grid = np.reshape(p_grid, (-1,params["extent"])) # convert back to 2d



##########

#extent = 100
#n_years = 1000
#
## for no spatial or temporal variance, set V = 0
#params_values = {"extent" : [extent], "n_years" : [n_years], 
#            
#            "temporal_G" : [0], # 0, 0.5, 1, 2, 5
#            "V_t" : [0],   "L_t" : [1], # 0, 1000, 2000, 3000, 4000    # 0.01, 0.1, 1, 10, 100, 1000
#            "fft_kernel_1D" : [RBF],
#
#            "spatial_G" : [0],  # 0, 10, 20, 30, 40
#            "V_s" : [0],  "L_s" : [1], # 0, 1000, 2000, 3000, 4000      # 0.01, 0.1, 1, 10, 100, 1000
#            "fft_kernel_2D" : [RBF],  
#
#            "lambda" : [10], "dispersal_frac" : [0.1], # 5, 10, 20, 40  # 0.1, 1, 5, 10
#            "density_dep" : [1], 
#            "T0" : [0], "GD" : [100], "Rmax" : [4]} # 7
#
## create df of all param combos
#params_df = params_dict2df(params_values)
#params_df = calc_params(params_df)
#
#params = params_df.loc[0].to_dict()
#
#t_grid = temp_spatial(params)
#pop = initial_population(t_grid, params)
#pop
#np.round(pop, 3)
#np.mean(pop,1)
#
#
#params["N"]
#
#R = pop_growth_rate(t_grid, params["T0"], params["GD"], params["Rmax"])
#R
#params["Rmax"]
#(R - 1) / (params["density_dep"] * R)
#initial_population(t_grid, params)
#params["N"]
#(params["Rmax"] - 1) / (params["density_dep"] * params["Rmax"]) 
#Rmax = 4
#(Rmax - 1) / (0.1 * Rmax) 
#
#B = (params["GD"]) / np.log(params["Rmax"]) # GD steepness, saturates at exp(GD / B)
#PGR = t_grid - params["T0"]
#PGR_alt = PGR / (1 + B * PGR)
#PGR[t_grid > params["T0"]] = PGR_alt[t_grid > params["T0"]] # temps > T0
#R =  np.exp(params["GD"] * PGR) 
#
#R_ = params["Rmax"] - 0.01
#dt = np.exp((params["GD"] * np.log(params["Rmax"])) / (1 - params["GD"] * np.log(params["Rmax"]) * B))
#
#
#
#t_grid = dt + params["T0"]
#
#

#1/(1 - params["N"] * params["density_dep"])


