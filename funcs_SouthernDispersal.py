#!/usr/bin/env python

###############################################################################################
# southern immigration
###############################################################################################

"""
Functions for calculation of dispersal from the southern border.
Done once at the start of each simulation 
"""

# imports 
import numpy as np
#import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

# functions 

def disp_kernel(x, y, lmbda): 
    """ Exponential dispersal kernel 
    Might want to be a bit more careful about normalisation 
        when on a lattice - the lambda^2 only normalises it to unity when
        integrated over 2D continuous space """
    #lmbda = params["lambda"]
    dist = np.sqrt(x**2 + y**2) 
    return lmbda**2 * np.exp(-lmbda * dist) / (2 * np.pi)


def plot_lmbda(lmbda=0.2):
    """lmbda the dispersal decay constant. larger lmbda, shorter dispersal distances"""
    dist = np.array(range(0,101))
    y = lmbda**2 * np.exp(-lmbda * dist) / (2 * np.pi)
    plt.plot(dist,y); plt.show()

def sum_kernel(x, y, dk, params):
    """ Works for a general dispersal kernel called by dk(x, y, ...)"""
    return sum(dk(x, y, params["lambda"]))

def generate_coordinate_ring(lout, lin):
    """ Computes a list of x, y coordinates for points within the region
    -lout <= x <= lout, -lout <= y <= -1, excluding the range 
    -lin <= x <= lin, -lin <= y <= -1
    Exception: if lin=1, then no inner range is excluded """
    x, y = [], []
    if lin > 0:
        for yyy in range(-1, -lin-1, -1):
            xx = list(range(-lout, -lin)) + list(range(lin+1, lout+1))
            yy = [yyy] * len(xx) # repeat yyy len(xx) times
            x = x + xx
            y = y + yy
        for yyy in range(-lin-1, -lout-1, -1): 
            xx = list(range(-lout, lout+1))
            yy = [yyy] * len(xx) 
            x = x + xx
            y = y + yy
    else:
        for yyy in range(-1, -lout-1, -1): 
            xx = list(range(-lout, lout+1))
            yy = [yyy] * len(xx) 
            x = x + xx
            y = y + yy
    return np.array(x), np.array(y)

def generate_coordinate_bracket(lout, lin, Y = -1):
    """ Generates a list of coordinates for a line between -lout <= x <=lout, 
    excluding -lin <= x <= lin
    Exception: if lin=0, then nothing excluded """
    #lout = 10 
    #lin = 1
    if lin > 0:
        x = list(range(-lout, -lin)) + list(range(lin+1, lout+1)) 
    else:
        x = list(range(-lout, lout+1))
    y = [Y] * len(x)
    return np.array(x), np.array(y)


def sum_kernel_gen(y_target, dk, params, coord, sk_previous = 0):
    """ general function that can sum kernel directly or recursively
    For direct: 
    coord =  generate_coordinate_ring
    sk_previous = 0 at start
    For recursive: 
        recursive = generate_coordinate_bracket """

    #y_target = params["extent"]
    #dk  = disp_kernel
    #coord = generate_coordinate_ring
    
    lout = 10
    tol = 1e-10 # relative tolerance
    x, y = coord(lout, 0) 
    sk_o = np.copy(sk_previous)
    sk_n = sk_o + sum_kernel(x, y-y_target,  dk, params) 
    #while (abs((sk_o - sk_n)/(sk_n)) > tol): # run time warning
    while (abs((sk_o - sk_n)) > tol * abs(sk_n)): # run time warning
        sk_o = np.copy(sk_n)
        lin = np.copy(lout)
        lout = 2 * lout
        #print(lout)
        x, y = coord(lout, lin)
        sk_n = sk_o + sum_kernel(x, y - y_target,  dk, params)
    return sk_n

def s_dispersal(params):
    """Dispersal potential into simulated landscape from the virtual landscape beyond the southern border"""
    sk = np.array(sum_kernel_gen(params["extent_y"]-1, disp_kernel, params, generate_coordinate_ring)) 
    sk_old = np.copy(sk)
    for Y in range(params["extent_y"]-2, -1, -1):
        # Y is distance from southern border of patch i in virtual landscape
        sk_new = sum_kernel_gen(Y, disp_kernel, params, generate_coordinate_bracket, sk_old) 
        sk = np.append(sk, sk_new)
        sk_old = np.copy(sk_new) 
    return sk


# need to specify disp_kernel as input 
def southern_immigration(params):
    """Immigration into simulated landscape from the virtual landscape beyond the southern border """
    sk = s_dispersal(params)
    sk[np.abs(sk) < params["tol"]] = 0 # set values lower than tol to zero 
    sk = np.repeat(sk[None], params["extent_x"], axis=1)
    #sk = sk.reshape((-1,params["extent"]))
    sk = np.flip(sk.reshape((-1,params["extent_x"])))
    # equilibrium pop density south of border * fraction of pop that disperse
    #imms = params["N"] * params["dispersal_frac"] * sk 
    return sk #imms 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def southern_immigration_new(kernel, params):
    """new method does not require the above code """
    wrap_size = params["wrap_size"]
    # create simulated landscape of zeros
    grid = np.zeros((params["extent_y"], params["extent_x"]))
    # add east/west zero padding
    padded_grid = np.pad(grid,((0,0),(wrap_size, wrap_size)), mode='constant')
    # add virtual landscape as padding (N * dispersal_frac) 
    padded_grid = np.pad(padded_grid,((wrap_size,0),(0,0)), mode='constant', constant_values = (params["N"] * params["dispersal_frac"])) 
    # convolve and cut off padding
    s_imm = fftconvolve(padded_grid, kernel, mode ="same")[wrap_size:,wrap_size:-wrap_size]
    s_imm[np.abs(s_imm) < params["tol"]] = 0 
    return s_imm[:,1]

