#!/usr/bin/env python

###############################################################################################
# temporal and spatial variation in temperature
###############################################################################################


"""
Functions for calculation smooth spatial and temporal temperture gradients 
and for spatial and temporal variation which involve stationary gaussian processes

"""
# imports
import numpy as np
from itertools import product
#import matplotlib.pyplot as plt

# functions 

def RBF(Q, L):
    """ squared exponential covariance """
    K = np.exp(-(Q)/(2*L**2))
    return K

def RBF_2(Q, L, V):
    """ squared exponential covariance """
    K = V * np.exp(-(Q)/(2*L**2))
    return K

def spatial_variation_fft(extent_x, extent_y, V, L, func):
    """ """
    #x1 = np.linspace(0, extent+1, extent+1, endpoint= False)
    #x2 = np.linspace(0, extent+1, extent+1, endpoint= False) 
    #x1x2 = np.array(list(product(x1, x2)))
    #D = np.sum(x1x2**2, axis=1).reshape((extent+1, extent+1))

    ##
    x1 = np.linspace(0, extent_y+1, extent_y+1, endpoint= False)
    x2 = np.linspace(0, extent_x+1, extent_x+1, endpoint= False)
    x1x2 = np.array(list(product(x1, x2)))
    # get distance from 0,0
    D = np.sum(x1x2**2, axis=1).reshape((extent_y+1, extent_x+1))
    ##

    # periodic boundaries
    Q4 = np.copy(D[:-1,:-1]) # bottom right quadrant
    Q3 = np.copy(np.flip(D,1)[:-1,:-1]) # bottom left
    Q2 = np.copy(np.flip(D)[:-1,:-1]) # top left
    Q1 = np.copy(np.flip(D,0)[:-1,:-1]) # top right

    # combined covariance functions
    K = func(Q4, L) + func(Q3, L) + func(Q2, L) + func(Q1, L)

    b = np.random.normal(size = (extent_y, extent_x))
    b2 = np.random.normal(size = (extent_y, extent_x)) * 1j
    Z = np.fft.ifft2((b + b2) * np.sqrt(np.abs(np.fft.fft2(K))))
    Z = V * Z.real
    return Z

#extent_y = 2**13 #- 552
#extent_x = 2**10
#wrap_size = 552 # can tolerate a lambda as low as 0.05
#spatial_temp_range = 10
#spatial_G = spatial_temp_range/(extent_y - wrap_size *3)
#
#np.random.seed(3); 
#zz = spatial_variation_fft(extent_x, extent_y, V = 4344, L = 100, func = RBF)
#np.min(zz); np.max(zz)
#plot_gp_2D(zz, extent_y, extent_x)
#
#t_grid = data["temp_grid"]
#p_grid = data["pop_grid"]
#results["Pop"]["200"]
#
#
#np.mean(zz)
#np.sqrt(np.var(zz))
#
#t_grid_bin = t_grid > 0 
#int(t_grid_bin)
#
#
#plt.imshow(t_grid, cmap='RdBu_r')
#plt.imshow(t_grid_bin, cmap='RdBu_r')
#plt.imshow(p_grid, cmap='RdBu_r')
#plt.colorbar()
#plt.show()
#
#plt.gca().invert_yaxis()
#plt.xlabel('Longitude ($x$)', fontsize=sx); plt.ylabel('Latitude ($y$)', fontsize=sy)
#plt.tick_params(axis='both', which='major', labelsize=12)
#plt.ylim(0,extent_y); plt.xlim(0,extent_x)
#plt.title("Spatial Gradient", loc='left', fontsize=sh,fontweight='bold')
#


def temporal_variation_fft(n_years, V, L, func):
    """   """
    x1 = np.linspace(0, n_years -1, n_years)
    x2 = np.flip(np.linspace(1, n_years, n_years))

    K = func(x1**2, L) + func(x2**2, L)
   
    b = np.random.normal(size = n_years)
    b2 = np.random.normal(size = n_years) * 1j 
    Z = np.fft.ifft((b + b2) * np.sqrt(np.abs(np.fft.fft(K))))
    Z = V * Z.real 
    return Z

def temporal_variation_fft_2(n_years, V, L, func):
    """   """
    x1 = np.linspace(0, n_years -1, n_years)
    x2 = np.flip(np.linspace(1, n_years, n_years))

    K = func(x1**2, L, V) + func(x2**2, L, V)
   
    b = np.random.normal(size = n_years)
    b2 = np.random.normal(size = n_years) * 1j 
    Z = np.fft.ifft((b + b2) * np.sqrt(np.abs(np.fft.fft(K))))
    Z = Z.real 
    return Z

# ~~~ does changing the length scale change the variance 

#import numpy as np
#np.random.seed(np.random.randint(10)); 
#
#
#A = 1.5 * np.sqrt(3000)
#
#vars = []
#for i in range(10):
#    vars += [np.mean(temporal_variation_fft(n_years=3000, V = A, L = 10, func = RBF)**2)]
#np.mean(vars)
#
#A = 1.5 * np.sqrt(extent_x * extent_y)
#
#extent_y = 2**13 #- 552
#extent_x = 2**10
#
#vars = []
#for i in range(10):
#    vars += [np.mean(spatial_variation_fft(extent_x, extent_y, V = A, L = 10, func = RBF)**2)]
#np.mean(vars)
#
#plt.plot(range(0,10000), temporal_variation_fft(n_years=10000, V = 100, L = 10, func = RBF), color ="red");plt.show()
#
#plt.hist(temporal_variation_fft(n_years=10000, V = 100, L = 10, func = RBF));plt.show()
#
## create df to store results 
#import pandas as pd
#import matplotlib.pyplot as plt
#df = pd.DataFrame(columns=['func','L','V','mean_var', 'mean_sd','n_years'])
#
#counter=0
#for funct in ["RBF", "RBF_2"]: #, "RBF_2"
#    for n_years in [10000, 50000]: #, 20000
#        for A in [10, 100, 200]: 
#            for L in [1, 10, 100, 1000]: #
#                vars = np.array([])
#                for i in range(1000):
#                    if funct == "RBF":
#                        v = np.mean(temporal_variation_fft(n_years=n_years, V = A, L = L, func = RBF)**2)
#                    elif funct == "RBF_2":
#                        v = np.mean(temporal_variation_fft_2(n_years=n_years, V = A, L = L, func = RBF_2)**2)
#                    vars = np.append(vars, v)
#                df.loc[counter,"func"] = funct
#                df.loc[counter,"n_years"] = n_years    
#                df.loc[counter,"V"] = A
#                df.loc[counter,"L"] = L
#                df.loc[counter,"mean_var"] = np.mean(vars)
#                df.loc[counter,"mean_sd"] = np.mean(np.sqrt(vars))
#                counter += 1
#                print(counter)
#     
#
#plt.plot(range(0,100000), temporal_variation_fft(n_years=100000, V = 10, L = 100, func = RBF), color ="red");plt.show()
#
#df.to_csv("../Results/spatio_temp_params_study2.csv")
#
#
#a = np.array([])
#for i in range(100):
#        b=np.sqrt(np.var(temporal_variation_fft(n_years=10000, V = 100, L = 100, func = RBF)))
#        a = np.append(a,b)
#np.mean(a)

# do this for temp and spatial 

#plt.plot(range(20000), (np.asarray(range(20000)) * 0.015) + temporal_variation_fft(n_years=20000, V = 100, L = 500, func = RBF));plt.show()

#np.random.seed(10)
#np.mean(temporal_variation_fft(n_years =20000, V = 300, L = 1000, func = RBF))

def temp_temporal(start, end, params):
    """ temporal change of temporature """
    temps = np.array(range(start, end)) * params["temporal_G"] 
    #temps = temps + temporal_variation_fft(params["n_years"], params["V_t"], params["L_t"], RBF)
    return temps


#plt.plot(range(0, 20000), np.array(range(0, 20000)) * temporal_G[3]);plt.show()


def temp_spatial(params):
    """ spatial grid of temperature """
    y = np.linspace(0, params["extent"], params["extent"], endpoint=False)
    yy = np.tile(y, (params["extent"],1)).T 
    t_grid = params["intercept"] - (yy * params["spatial_G"])
    # if temp == T0 across landscape and wont change with time, set to be > T0
    if params["spatial_G"] == 0:# or params["temporal_G"] == 0:
        t_grid = 20
    t_grid = t_grid + spatial_variation_fft(params["extent"], params["V_s"], params["L_s"], RBF)
    return t_grid


#extent_y = 2000
#extent_x = 1
#T0 = 0 
#spatial_G = 0.02
##
#y0 = extent_y - 200  # starting position of climate range margin
### ensures initial temp at y0 == T0
#intercept = T0 + (spatial_G * y0) 
##
#y = np.linspace(0, extent_y, extent_y, endpoint=False)
#lats = np.tile(y, (extent_x, 1)).T 
#Tt = intercept - (lats * spatial_G)
#Tt[:,0]

def temp_spatial_lats(extent_y):
    """ latitude grid"""
    y = np.asarray(range(extent_y))
    #yy = np.tile(y, (extent_x, 1)).T 
    return y#y


def temp_spatial_grad(params, lats):    
    """ take lats and create spatial gradient of temp"""
    #lats = data['lat_grid']
    #print(lats)
    t = params["intercept"] - (lats * params["spatial_G"])
    t_grid = np.tile(t, (params["extent_x"],1)).T

    # if temp == T0 across landscape and wont change with time, set to be > T0
    if params["spatial_G"] == 0:# or params["temporal_G"] == 0:
        t_grid = 1000 # high enough to ensure max growth rate
    return t_grid 

#params["spatial_G"] = 30
#params["intercept"] = params["T0"] + (params["spatial_G"] * params["y0"]) 
#(params["intercept"] - (lats * params["spatial_G"]))[int(params["y0"])]
#(params["intercept"] - (lats * params["spatial_G"]))



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rolling landscape functions

# update it 
def update_pop_grid(pop_grid, n):
    """ """
    #n = 1
    # remove rows
    #a = data['pop_grid'][n:,]
    ## create new rows
    #row = np.zeros(params['extent_x'])
    #rows = np.repeat(row, n).reshape(n, params['extent_x'])
    ## update original 
    #data['pop_grid'] = np.append(a,rows, 0)

    pop_grid[:-n] = pop_grid[n:]
    pop_grid[-n:] = 0
    #return data['pop_grid']


def update_lat_grid(lat_grid, n):
    """ """
    # remove rows
    #a = data['lat_grid'][n:,]
    # create new rows
    #next_lat = int(np.max(a)) + 1
    #lats = np.asarray(range(next_lat, next_lat + n))
    #lat_grid = np.tile(lats.astype(float), (params['extent_x'],1)).T
    ## update original 
    #data['lat_grid'] = np.append(a,lat_grid, 0)

    lat_grid[:-n] = lat_grid[n:]
    lat_grid[-n:] += n
    #return data['lat_grid'] 


def update_temp_var(temp_var_grid, n):
    """ """
    #n = 10
    # Store rows
    #rows = data['temp_var_grid'][:n,:]
    ## remove rows
    #a = data['temp_var_grid'][n:,]
    ## swap rows
    #data['temp_var_grid'] = np.append(a, rows, 0)

    rows = np.copy(temp_var_grid[:n, :])
    temp_var_grid[:-n, :] = temp_var_grid[n:, :]
    temp_var_grid[-n:, :] = rows
    #return data['temp_var_grid']

#plot_gp_2D(data['temp_var_grid'], 100,100)


####### Plots

def plot_gp_2D(Zp, extent_x, extent_y):
    """ only works if one prior sampled """
    #x1 = np.linspace(0, extent, extent, endpoint= False)
    #x2 = np.linspace(0, extent, extent, endpoint= False) 
    x1 = np.linspace(0, extent_y, extent_y, endpoint= False)
    x2 = np.linspace(0, extent_x, extent_x, endpoint= False)

    x1x2 = np.array(list(product(x1, x2)))
    #X0p, X1p = x1x2[:,0].reshape((extent,extent)), x1x2[:,1].reshape((extent,extent)) # plotting space
    X0p, X1p = x1x2[:,0].reshape((extent_y,extent_x)), x1x2[:,1].reshape((extent_y,extent_x)) # plotting space
    fig = plt.figure(figsize=(10,8)) # surface plot 
    ax=fig.add_subplot(111)
    plot = ax.pcolormesh(X0p, X1p, Zp) # contour plot
    fig.colorbar(plot)
    #ax = fig.add_subplot(111, projection='3d')   # surface plot 
    #surf = ax.plot_surface(X0p, X1p, Zp, rstride=1, cstride=1, cmap='jet', linewidth=0, antialiased=False) # surface plot
    plt.show()

#plot_gp_2D(samples_2D, x1x2, extent)
