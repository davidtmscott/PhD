#!/usr/bin/env python

###############################################################################################
# collect results
###############################################################################################

"""
Functions to collect results at each iteration of the simulation

Results include the population range margin and the climatic range margin 

"""

# imports
import numpy as np
import pandas as pd 

# import functions 
from funcs_PopDynamics import carrying_capacity

# functions 

def create_results_dic(params):
    """ return dict of dicts to store results for one parameter set"""
    # limits, temps, landcover of length n_years
    # R and pop of length extent_y
    results = { "limits" : {"P_lim_k" : np.zeros(params["n_years"], dtype="float64"), 
                            #"P_lim_D0.01" : np.zeros(params["n_years"], dtype="float64"), 
                            "P_lim_D" : np.zeros(params["n_years"], dtype="float64"), 
                            "C_lim_lin": np.zeros(params["n_years"], dtype="float64")
                            }, 
                "temps" : [], 
                "landcover" : {"can_gen_pop" : np.zeros(params["n_years"], dtype="float64"), 
                               "occupied" : np.zeros(params["n_years"], dtype="float64")}, 
                "R"  : dict(), "Pop" : dict(),
                "lats" : dict(), "K" : dict()} # 
    return results


def interpolation_density(params, p_grid, l_grid, threshold): 
    """ calculates exact position of species range margin 
    using a density threshold """
    #l_grid = data['lat_grid']
    #p_grid = data['pop_grid']
    #threshold= 0.01
    m_grid = np.mean(p_grid, 1)
    try:
        #threshold=0.01
        #y0 = params["N"] * threshold
        y0 = threshold
        # where pop > threshold
        #x1_ = np.where(m_grid > y0)[0][-1] # had this when y0 was a frac of N
        x1_ = np.where(m_grid >= y0)[0][-1] 
    except:
        x0 = float("NaN")
    else:
        x2 = x1_ + 1        
        #print(x1_)
        y1 = m_grid[x1_]
        #print(y1) 
        y2 = m_grid[x2]
        # get actual latitude (for shifting landscape)
        x1 = l_grid[x1_] #,0 
        x2 = x1 + 1
        # exponential 
        k = (np.log(y2) - np.log(y1)) / (x1 - x2)
        a = k * x1 + np.log(y1)
        x0 = (-np.log(y0) + a) / k
        # linear
        #b = (y2 - y1) / (x2 - x1)
        #a = y1 - (b * x1)
        #x0 = (y0 - a) / b

        #A = np.array([[1, x1],[1, x2]])
        #B = np.array([y1, y2])
        #np.linalg.solve(A,B)
    return x0 

#interpolation_density(params, data["pop_grid"], data["lat_grid"], 0.01)

def integrate_K(K, total_N):
    """ find latitude where K == N """
    lat = 0
    total_K = 0
    while (total_K < total_N):    
        # sum total K 
        lat_K = np.sum(K[lat])
        total_K = np.add(total_K, lat_K)
        lat += 1
    return lat - 1


def interpolate_K(total_N, K, K_lat, l_grid):
    """ finds exact location of bio range margin based on 
    K carry capacity of landscape"""
    y0 = abs(total_N)
    x1 = K_lat - 1
    y1 = np.sum(K[:x1+1])    
    x2 = K_lat
    y2 = np.sum(K[:x2+1])

    # get actual latitude (for shifting landscape)
    x2 = l_grid[K_lat] #,0 
    x1 = x2 - 1
    
    # exponential
    k = (np.log(y2) - np.log(y1)) / (x1 - x2)
    a = k * x1 + np.log(y1)
    x0 = (-np.log(y0) + a) / k

    # linear
    #b = (y2 - y1) / (x2 - x1)
    #a = y1 - (b * x1)
    #x0 = (y0 - a) / b
    return x0 



def range_margin_K(params, data):
    """ appends nothern latitude of range margin to results_dict """
    #p_grid = data['pop_grid']   
    K = carrying_capacity(data["R"], params)
    total_N = np.sum(data['pop_grid'])
    #print(data['pop_grid'])
    # breaks here in the case where total_N > total_K
    #if total_N <= np.sum(K):
    try:
        K_lat = integrate_K(K, total_N) 
        range_lat = interpolate_K(total_N, K, K_lat, data["lat_grid"]) #+ lat_counter 
    except:
        range_lat = np.nan
    return range_lat


def climate_margin_linear(params, temp):
    """ latitude of climate margin. where temp = T0 at time point j"""  
    y_lat = (params["intercept"] - params["T0"] + temp) / params["spatial_G"]
    return y_lat

def climate_margin_T0(params, data, j):
    """ latitude of climate margin. where temp = T0 at time point j"""  
    T = data['temp_grid'] + data['temps'][j] # update temperature   
    lat_temps = np.mean(T, 1)
    # position of last temp > T0
    y_lat = int(np.where(lat_temps >= params["T0"])[0][-1:]) 
    #if lat_temps[y_lat] > 0:
    # interpolate
    y0 = params["T0"]
    x1 = y_lat
    y1 = lat_temps[x1]
    x2 = y_lat + 1
    y2 = lat_temps[x2]
    # get actual latitude (for shifting landscape)
    x1 = data["lat_grid"][y_lat] #,0 
    x2 = x1 + 1
    
    # exponential
    #k = (np.log(y2) - np.log(y1)) / (x1 - x2)
    #a = k * x1 + np.log(y1)
    #y_lat = (-np.log(y0) + a) / k

    # linear
    b = (y2 - y1) / (x2 - x1)
    a = y1 - (b * x1)
    y_lat = (y0 - a) / b
    return y_lat


def climate_margins(params, data, dic, j):
    """ latitude of climate margin. where temp = T0 at time point j"""  
    if params["spatial_G"] != 0:
        dic["limits"]["C_lim_lin"][j] = climate_margin_linear(params, data['temps'][j])
        #try:
        #    dic["limits"]["C_lim_T0"][j] = climate_margin_T0(params, data, j)
        #except:
        #    dic["limits"]["C_lim_T0"][j] = float("NaN")    
    else: # no climate margin if spatial_G is 0
        dic["limits"]["C_lim_lin"][j] = float("NaN")
        #dic["limits"]["C_lim_T0"][j] = float("NaN")
    #return dic

# store less densities 
# dont store K 
# only collect results after burn in e.g. 1000 time steps 

def collect_results(dic, params, data, R, j):
    """ """
    #dic = range_margin_density(dic, params, data['pop_grid'], data["lat_grid"]) 
    #dic["limits"]["P_lim_k"][j] = range_margin_K(params, data)
    #dic["limits"]["P_lim_D0.01"][j] = interpolation_density(params, data["pop_grid"], data["lat_grid"], 0.01)
    dic["limits"]["P_lim_D"][j] = interpolation_density(params, data["pop_grid"], data["lat_grid"], 0.014352285636780126)


    # climate margins
    climate_margins(params, data, dic, j) #dic = 
    # % land cover habitat v non habitat - binary R>=1 is habitat
    #dic["landcover"]["can_gen_pop"][j] = np.sum((R >= 1).view(np.int8)) / (params["extent_y"] * params["extent_x"])
    # % land cover occupied
    #dic["landcover"]["occupied"][j] = np.sum((data["pop_grid"] > 0).view(np.int8)) / (params["extent_y"] * params["extent_x"])
    # take snapshot every 100 years
    if j == 0 or j % 1 == 0:
        dic["Pop"][str(j)] = np.mean(data['pop_grid'], axis = 1)  
        dic["lats"][str(j)] = data['lat_grid'].copy() #np.mean(data['lat_grid'], axis = 1)
        dic["R"][str(j)] = np.mean(R, axis = 1)
        dic["K"][str(j)] = np.mean(carrying_capacity(R, params), axis = 1)
        # ~~ could save data["lat_grid"] here
    #return dic

def results_dict2df(dic, params, x, col):
    """ convert final results dicts to dfs for saving """
    #dic=res_dic
    #dic = results["limits"]
    #x=n_years
    #col = "years"
    ps = pd.DataFrame(params, index = [0])
    df = pd.DataFrame.from_dict(dic)
    # repeat each combo of param values x times in a df
    #x=10
    params_results = pd.DataFrame(np.repeat(ps.values, x, axis = 0))
    params_results.columns = ps.columns
    # combine so each result is in a row with correspecing combo of param values
    res = pd.concat((params_results, df), axis = 1)
    # add in years or space variable as column
    res[str(col)] = np.tile(np.array(range(0, x)), ps.shape[0])
    return res


def format_results(params, results):
    margins = results_dict2df(results["limits"], params, params["n_years"], "year")
    pop_dens = results_dict2df(results["Pop"], params, params["extent_y"], "space") # params_df.iloc[[i]]
    Rs = results_dict2df(results["R"], params, params["extent_y"], "space")
    #Ks = results_dict2df(results["K"], params, params["extent_y"], "space")
    # calc lag 
    #margins["lag_k"] =  margins["P_lim_k"] - margins["C_lim"]
    #margins["lag_d"] =  margins["P_lim_d"] - margins["C_lim"]
    # calc invasion speed 
    p_lim = margins["P_lim_D"] #["P_lim_D0.01"]
    speed = (p_lim.iloc[-1] - p_lim.iloc[50]) / (len(p_lim) - 50)
    return margins, pop_dens, Rs, speed #Ks

#def pop_density(dic, p_grid, j):
#    """ mean population densities at certain time points """
#    if (j in list(map(int, dic.keys()))):
#        dic[str(j)].extend(np.mean(p_grid, axis = 1))
#    return dic

#def dens_dict(n_years):
#    """ """
#    a = np.linspace(0,n_years + 100, int(n_years / 100 + 1), endpoint=False, dtype=int)
#    #a = np.array(range(0,10))
#    a[-1] -= 1
#    d = {str(key) : [] for key in a}
#    return d


##################
# test interpolation 


#extent = 10
##np.random.seed(1)
#
#m_grid = np.random.exponential(2,10)
#m_grid[::-1].sort()
#
#space = list(range(0,extent))
#
#space_interpolated = np.array([])
#y0_array = np.array([])
#
#for x1 in space:
#    #x1 = 0
#    y1 = m_grid[x1]    
#    
#    x2 = x1 + 1
#    if x2 == extent:
#        space_interpolated = np.append(space_interpolated, float("NaN"))
#        y0_array = np.append(y0_array, float("NaN"))
#        break
#    y2 = m_grid[x2]
#    
#    y0 = np.random.uniform(y1,y2)
#    
#    #print("neww.......")
#    #k = -y2 / x2*x1*np.log(y1)
#    #a = y1/np.exp(-k*x1)
#    #a = y1 * np.exp(-y2/(x2 * np.log(y1)))
#    #print(k)
#    #print(a)
#    #x0 = (-np.log(y0) + np.log(a)) / k
#    
#    #print("old..........")
#    k = (np.log(y2) - np.log(y1)) / (x1 - x2)
#    #k = np.log(y2 - y1) / (x2 - x1)
#    a = np.exp(k * x1 + np.log(y1))
#    #print(k)
#    #print(a)
#    x0 = (-np.log(y0) + np.log(a)) / k
#    
#    #b = (y2 - y1) / (x2 - x1)
#    #a = y1 - (b * x1)
#    #x0 = (y0 - a) / b
#
#    print("jjjjj")
#    print(y0)
#    print(x0)
#    space_interpolated = np.append(space_interpolated, x0)
#    y0_array = np.append(y0_array, y0)
#
#import matplotlib.pyplot as plt
#
#
#plt.plot(space, m_grid)
#plt.scatter(space_interpolated, y0_array)
##plt.xlim(left = 0, right=10)
#plt.show()
#
#