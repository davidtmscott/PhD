#!/usr/bin/env python

###############################################################################################
# dispersal
###############################################################################################

"""
Functions to calculate dispersal 

includes dispersal kernel , disptance calculation and fft convultions
"""
# imports
import numpy as np
import pyfftw

#import matplotlib.pyplot as plt

# functions

def kernel_uniform(distance):
    """ distance = how many cells it can disperse over
    assumes array with odd number of elements
    centre cell = 0, all other cells sum to 1 """
    cells = (distance * 2) + 1
    kernel = np.random.random_integers(0, cells**2 ,(cells,cells))
    np.put(kernel, kernel.size // 2, 0) # middle = 0
    kernel = kernel / kernel.sum()
    return kernel


def distances(extent_x):
    """ I removed extent_y """
    size_x = extent_x * 2 - 1
    size_y = extent_x * 2 - 1
    x_arr, y_arr = np.mgrid[0:size_y, 0:size_x]  
    centre_x = np.floor(size_x/2)
    centre_y = np.floor(size_y/2)  
    dists = np.sqrt((x_arr - centre_y)**2 + (y_arr -  centre_x)**2)
    return dists

def kernel_neg_exp(extent_x, lmbda):
    """ negative exponential kernel. centre cell = 0, all other cells sum to 1
    kernel size = extent * 2 - 1 """
    #extent_x = 1000
    #extent_y = 1000
    #lmbda=0.1
    dists = distances(extent_x)
    kernel = lmbda**2 * np.exp(-lmbda * dists) / (2 * np.pi)
    np.put(kernel, kernel.size // 2, 0) # centre = 0
    kernel = kernel / kernel.sum() # convert to probs
    return kernel


def kernel_neg_exp_new(extent_x, lmbda):
    """ negative exponential kernel. centre cell = 0, all other cells sum to 1
    kernel size = extent * 2 - 1 """
    #extent_x = 10
    #extent_y = 10
    #lmbda=0.1
    dists = distances(extent_x)
    kernel = lmbda**2 * np.exp(-lmbda * dists) / (2 * np.pi)
    #np.put(kernel, kernel.size // 2, 0) # centre = 0
    #kernel = kernel / kernel.sum() # convert to probs
    
    # ~~ set tol cut off
    kernel[kernel < 1e-10]  = 0
    # ~~ make all non centre cells sum to 1 
    non_centre_cells = np.delete(kernel, kernel.size // 2) # select all but centre cell
    non_centre_cells = non_centre_cells / non_centre_cells.sum() # have them all sum to 1
    # insert centre as 0 by splitting in two halfs
    half1 = non_centre_cells[:non_centre_cells.size // 2]# first half 
    half2 = non_centre_cells[non_centre_cells.size // 2:]# second half 
    centre = np.asarray([0]) # 0 centre 
    kernel_new = np.concatenate((half1, centre, half2))  
    # resahpe to 2D kernel
    kernel_new = np.reshape(kernel_new, (-1, kernel.shape[1]))
    return kernel_new    


def cut_kernel(Kyx):
    """ remove rows and cols containing only zero """ 
    Kyx = Kyx[:, ~np.all(Kyx == 0, axis=0)]
    Kyx = Kyx[~np.all(Kyx == 0, axis=1), :]
    return Kyx

def fft_kernel(kernel, wrap_size, extent_x, extent_y):
    """ precompute the FFT of the kernel. this is only done once at start"""
    # pad the grid with periodic boundary conditions
    #padded_grid = np.pad(grid,((wrap_size,wrap_size),(0,0)), mode='constant')
    # Compute the size of the padding needed to make the filter the same size as the signal
    #padding_size = [(padded_grid.shape[i]-kernel.shape[i]) for i in range(2)]
    # ~
    #padding_size = [(extent_y + (wrap_size * 2)) - kernel.shape[0], (extent_x) - kernel.shape[0]]
    padding_size = [(extent_y) - kernel.shape[0], (extent_x) - kernel.shape[0]] # (extent_y + 552)
    # Compute the amount of padding needed on each side of the filter
    pad_width = [(int(np.ceil(padding_size[i]/2)), int(np.floor(padding_size[i]/2))) for i in range(2)]
    # Pad the filter with zeros using the pad_width
    kernel_padded = np.pad(kernel, pad_width, mode='constant', constant_values=0)
    # Shift the padded filter using the fftshift function to center it
    kernel_padded_centered = np.fft.fftshift(kernel_padded)
    kernel_fft = np.fft.rfft2(kernel_padded_centered)
    return kernel_fft



def make_kernel(params):
    """ I changed this to only use extent_x, 
    if its big enough for extent_x its big enough for extent_y as 
    x will either be equal to or less than y """
    Kernel = kernel_neg_exp_new(params["extent_x"], params["lambda"])
    Kernel = cut_kernel(Kernel)
    # get wrap_size 
    params["wrap_size"] = Kernel.shape[0] // 2
    # get fft of kernel
    kernel_fft = fft_kernel(Kernel, params["wrap_size"], params["extent_x"], params["extent_y"])
    return Kernel, kernel_fft, params



# plot kernel 
def plot_kernel(extent_x, extent_y, lmbda):
    """ function to plot 3D kernel"""
    X = np.linspace(-extent_x,extent_x,extent_x*2-1)
    Y = np.linspace(-extent_y+1,extent_y,extent_y*2-1)
    X, Y = np.meshgrid(X,Y)
    Z = kernel_neg_exp(extent_x, extent_y, lmbda)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(X, Y, Z, cmap="plasma")
    ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
    plt.savefig('../Poster/figures/dispersal.png')
    plt.show()

#plot_kernel(100, 100, 0.01)


def emigration(p_grid, dispersal_frac):
    """ can change this to calculate an individual dispersal prob for each cell or, 
    to calculate a random dispersal prob each time function is called """
    dispersers = p_grid * dispersal_frac
    return dispersers

from scipy.signal import fftconvolve

def immigration(dispersers, kernel, tol):
    """ distributes dispersers as immigrants based on dispersal kernel
    calls convolve from scipy.ndimage
    E/W immigration is periodic ("wrap")
    N/S emmigration is absorbing ("constant"), so its padded
    S immigration is calc from second function 
    pads input with zeros. 
    use oaconvolve from scipy.signal. It uses fft. Better than using fftconvolve as my kernel is twice the size of landscape
    no need for padding, just subset and wrap """    
    #kernel = data['kernel']
    #tol = params["tol"]
    wrap_size = kernel.shape[0] // 2
    # ~~ remove north and south 
    #full = oaconvolve(dispersers, kernel, mode='full', axes=None)[wrap_size:-wrap_size,:]
    full = fftconvolve(dispersers, kernel, mode='full', axes=None)[wrap_size:-wrap_size,:]
    full[np.abs(full) < tol] = 0 # set cut off tol 1e-15
    # ~~ select immigrants within landscape 
    same = full[:, wrap_size:-wrap_size]
    # ~~ add to the east indivds who dispersed west 
    west_dispersers = full[:, 0: wrap_size]
    same[:,-west_dispersers.shape[1]:] = same[:,-west_dispersers.shape[1]:] + west_dispersers
    # ~~ add to the west individs who disperser east
    east_dispersers = full[:, -wrap_size:]  
    same[:,:east_dispersers.shape[1]] = same[:,:east_dispersers.shape[1]] + east_dispersers  
    # ~~~ alternative 
    #padded_grid = np.pad(dispersers, ((0,0), (wrap_size,wrap_size)), mode= "wrap")
    #same = fftconvolve(padded_grid, kernel, mode='same')[:,wrap_size:-wrap_size]
    #same[np.abs(same) < tol] = 0 
    return same  


def immigration_new(grid, kernel_fft, wrap_size, tol):
    #padded_grid = np.pad(grid,((wrap_size,wrap_size),(0,0)), mode='constant')
    padded_grid = np.pad(grid,((552,0),(0,0)), mode='constant')
    # compute the fft for pop grid
    #padded_grid_fft = pyfftw.interfaces.numpy_fft.rfft2(padded_grid)
    # apply the convolution using FFT-based method
    #convolved = pyfftw.interfaces.numpy_fft.irfft2(padded_grid_fft * kernel_fft)[wrap_size:-wrap_size,:]
    convolved = pyfftw.interfaces.numpy_fft.irfft2(pyfftw.interfaces.numpy_fft.rfft2(padded_grid) * kernel_fft)[552:,:]
    convolved[np.abs(convolved) < tol] = 0 # set cut off
    return convolved 


def immigration_new2(grid, kernel_fft, tol, dispersal_frac):
    padded_grid = np.pad(grid * dispersal_frac,((552,0),(0,0)), mode='constant')
    # apply the convolution using FFT-based method
    convolved = pyfftw.interfaces.numpy_fft.irfft2(pyfftw.interfaces.numpy_fft.rfft2(padded_grid) * kernel_fft)[552:,:]
    convolved[np.abs(convolved) < tol] = 0 # set cut off
    np.multiply(grid, (1 - dispersal_frac), out=grid)
    grid += convolved
    return grid

def immigration_new2_class_object(grid, kernel_fft, tol, dispersal_frac, data):
    non_disperser = grid * (1 - dispersal_frac)
    grid *= dispersal_frac
    # apply the convolution using FFT-based method
    b = data["fft"]()  #fft_object()
    b *= kernel_fft
    grid = data["ifft"]()  #ifft_object()
    grid[np.abs(grid) < tol] = 0 # set cut off
    grid += non_disperser
    #return grid


import threading

#@profile
def dispersal(params, data):
    """ first selct the number of disersers in emigration matrix 
    subtract emigration matrix (dispersers) from pop grid 
    convolve with dispersal kernel to new immigration matrix 
    add immigration matrix to pop_grid """
    #p_grid = data['pop_grid']
    #dispersal_frac = params['dispersal_frac']
    #kernel = data["kernel"]
    
    #dispersers = emigration(p_grid, params["dispersal_frac"])
    #p_grid -= dispersers
    #p_grid += immigration_new(dispersers, kernel_fft, params["wrap_size"], params["tol"])

    #p_grid = immigration_new2(p_grid, kernel_fft, params["tol"], params["dispersal_frac"])
    #immigration_new2_class_object(grid, kernel_fft, params["tol"], params["dispersal_frac"], data)

    non_disperser = data["pop_grid"] * (1 - params["dispersal_frac"])
    data["pop_grid"] *= params["dispersal_frac"]
    # apply the convolution using FFT-based method
    b = data["fft"]()  #fft_object()
    b *= data["kernel_fft"]
    data["pop_grid"] = data["ifft"]()  #ifft_object()
    data["pop_grid"][np.abs(data["pop_grid"]) < params["tol"]] = 0 # set cut off
    data["pop_grid"][-200:,:] = 0
    data["pop_grid"] += non_disperser
    data["pop_grid"] += data["southern_immigrants"][:, np.newaxis]
    #return grid

#   10    2.175    0.217   62.114    6.211 M:\PhD\project2\Code\funcs_Dispersal.py:182(dispersal)
# 159.449 seconds


# test and see if dispersal changes pop


# code to inspect southern immigration
# get length of kernel 
# create 2d lattive length of dispersal kernel 
# pop at k 
# attach to pop grid  
# set pop grid to zero 
# get dispersers

#
#p_grid = data['pop_grid']
#dispersal_frac = params['dispersal_frac']
#kernel = data["kernel"]
#tol = params["tol"]
#dispersers = emigration(p_grid, params["dispersal_frac"])
#
#
#wrap_size = kernel.shape[0] // 2
#
#virtual_pop = np.full((wrap_size,1000),(params["N"] * dispersal_frac))
#
#p_grid.fill(0) 
#
#
#pop = np.concatenate((virtual_pop,p_grid), axis=0)
#
#same = fftconvolve(pop, kernel, mode='same')[wrap_size:,:]
#same[np.abs(same) < tol] = 0 
#
#same.shape
#params
#same[1,500]
#
#southern_immigrants[1,0]
#
## ~~~
## code to check density of pop at southern border. if at equilibrium.
#(params["N"] - p_grid[0,0])/params["N"]
#p_grid[0,0]/params["N"]
#plt.plot(range(1000), p_grid[:,0]);plt.show()
#
#
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this works 
## ~~~~~~~~~~~~~~~~~~~
#
#import time 
#
#extent_y = 1000
#extent_x = 1000
#
#p_grid = data['pop_grid']
#dispersal_frac = params['dispersal_frac']
#kernel = data["kernel"]
#tol = params["tol"]
#grid = emigration(p_grid, params["dispersal_frac"])
#
#convolved[150:250,100]
#cc[250:350,100]
#grid[150:250,100]
#
## create a 2D lattice grid and a 2D dispersal kernel
#y = np.linspace(0, extent_y, extent_y, endpoint=False)
#grid = np.tile(y, (extent_x, 1)).T 
#
#
#grid = np.zeros((10, 10))
##grid = np.ones((10, 10))
#
#grid[1,0] = 1
##grid.fill(1)
##grid[:] = list(range(10))
##kernel = np.ones((5, 5)) / 25
#wrap_size = kernel.shape[0]//2
#
#kernel[wrap_size,wrap_size] = 0
#
#
#extent_y = 1000
#extent_x = 1000
#
#p_grid = data['pop_grid']
#dispersal_frac = params['dispersal_frac']
#kernel = data["kernel"]
#tol = params["tol"]
#grid = emigration(p_grid, params["dispersal_frac"])
#
#wrap_size = kernel.shape[0]//2
#
## perform fft on kernel once, then save
#kernel_fft = fft_kernel(grid,kernel,wrap_size)
#
## calculate immigration from grid of dispersers
#dispersers = immigration_new(grid, wrap_size, kernel_fft)
#
#
#ss = southern_immigration_new(kernel, params, wrap_size)
#
#
#dispersers + ss
#
#

#
#print(np.allclose(convolved, convolved_grid))
#
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## southern immigration
#
## create grid all zeros
## create kernel
## east west pad grid with zeros length of wrap_size
## create virtual grid 
## same ncols as grid 
## nrows as wrap_size
## fill to carrying capacity 
## multiply by dispersal fraction
## perform convolution 
## discard virtual landscape, and padding 
#
#extent_y = 10000
#extent_x = 1000
#
#p_grid = data['pop_grid']
#dispersal_frac = params['dispersal_frac']
#kernel = data["kernel"]
#tol = params["tol"]
#grid = emigration(p_grid, params["dispersal_frac"])
#
#
#grid = np.zeros((10000, 1000))
#
#kernel = np.ones((200, 200)) / 25
#wrap_size = kernel.shape[0]//2
#kernel[wrap_size,wrap_size] = 0
#
#
#padded_grid = np.pad(grid,((0,0),(wrap_size,wrap_size)), mode='constant')
## add virtual landscape padding 
#padded_grid = np.pad(padded_grid,((wrap_size,0),(0,0)), mode='constant', constant_values = (params["N"] * params["dispersal_frac"])) 
##padded_grid = np.pad(padded_grid,((wrap_size,0),(0,0)), mode='constant', constant_values=1) * 0.2
#
#
#cc = fftconvolve(padded_grid, kernel, mode ="same")[wrap_size:,wrap_size:-wrap_size]
#cc[np.abs(cc) < tol] = 0 
#cc.shape
#cc[1,0]
#
#cc + convolved
#
#



#######################################

#def d_k(dists, lmbda): 
#    #dists=Q4
#    d = lmbda**2 * np.exp(-lmbda * dists) / (2 * np.pi) 
#    return d 
#
#extent = 10
#x = np.linspace(0, extent+1, extent+1, endpoint= False)
#y = np.linspace(0, extent+1, extent+1, endpoint= False) 
#xy = np.array(list(product(x, y)))
## get distance from 0,0
#D = np.sqrt(np.sum(xy**2, axis=1)).reshape((extent+1, extent+1))
##D = np.round(D, 2)
## periodic boundaries
#Q4 = np.copy(D[:-1,:-1]) # bottom right quadrant
#Q3 = np.copy(np.flip(D,1)[:-1,:-1]) # bottom left
#Q2 = np.copy(np.flip(D)[:-1,:-1]) # top left
#Q1 = np.copy(np.flip(D,0)[:-1,:-1]) # top right
#
#lmbda = 0.1
## combined covariance functions
#K = d_k(Q4, lmbda) + d_k(Q3, lmbda) + d_k(Q2, lmbda) + d_k(Q1, lmbda)
#b = np.random.normal(size = (extent, extent))
#b2 = np.random.normal(size = (extent, extent)) * 1j
#
#Z = np.fft.ifft2((b + b2) * np.sqrt(np.abs(np.fft.fft2(K))))
#Z = V * Z.real
#Z



######
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#extent = 10
#lmbda = 0.1
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#size = extent * 2 - 1
#x_arr, y_arr = np.mgrid[0:size, 0:size]  
#centre = np.floor(size/2)
#dists = np.sqrt((x_arr - centre)**2 + (y_arr -  centre)**2)
#
#kernel = lmbda**2 * np.exp(-lmbda * dists) / (2 * np.pi)
#np.put(kernel, kernel.size // 2, 0) # centre = 0
#
##kernel.ravel()[180]
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#p_grid = np.zeros((extent,extent))
#np.put(p_grid, (extent**2 // 2)+4, 1)
##(extent**2 // 2)
#
#dispersers = emigration(p_grid, 1)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
#wrap_size = int((kernel.shape[0] - 1) / 2)
#full = oaconvolve(dispersers, kernel, mode='full', axes=None)
#full = full[wrap_size:extent+wrap_size,:] # remove north and south 
#full[np.abs(full) < 1e-15] = 0 # set cut off tol
#same = full[:, wrap_size:extent+wrap_size] # select range 
#same[:,extent-wrap_size:] = same[:,extent-wrap_size:] + full[:,:wrap_size] # add to the east indivds who dispersed west   
#same[:,:wrap_size] = same[:,:wrap_size] + full[:, extent+wrap_size:] # add to west indivds that dispers to the east
#np.round(same, 8)
#####
#
#wrap_size = int((kernel.shape[0] - 1) / 2)
#full = oaconvolve(dispersers, kernel, mode='full', axes=None)[wrap_size:wrap_size*2+1,:] # remove north and south 
#full[np.abs(full) < 1e-15] = 0 # set cut off tol
#same = full[:, wrap_size:wrap_size*2+1] # select range 
#same[:,1:] = same[:,1:] + full[:, :wrap_size] # add to the east indivds who dispersed west   
#same[:,:-1] = same[:,:-1] + full[:, wrap_size*2+1:] 
##same
#np.round(same, 8)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#same[extent-wrap_size-1:,wrap_size+1:extent-2]
#same[wrap_size+1:extent-wrap_size,:wrap_size]
#
#
#wrap_size = int((kernel.shape[0] - 1) / 2)//2
#f = kernel[4:-5,:]
#s = f[:, 5:-4] # select range 
#s[:,-5:] = s[:,-5:] + f[:,:5] # add to the east indivds who dispersed west   
#s[:,:4] = s[:,:4] + f[:,-4:] 
#s
#

