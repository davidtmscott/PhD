#!/usr/bin/env python

"""
Code to calculate lmabdas for given dispersal lengths
"""

# imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from funcs_Dispersal import distances, kernel_neg_exp_new

# ~~~ find lambdas for specific dispersal distances
lmbda=0
positions = np.array([])
lambdas = np.array([])
for i in range(15000):
    lmbda += 0.0001
    # get y values 
    Kernel = kernel_neg_exp_new(2**10, lmbda)
    #Kernel = cut_kernel(Kernel)
    dists = distances(2**10)
    position = np.sum(Kernel * dists)
    print(position)
    print(i)
    positions = np.append(positions, position)
    lambdas = np.append(lambdas, lmbda)

plt.plot(lambdas, positions);plt.show()

df = pd.DataFrame(columns=['positions','lambdas'])

df["lambdas"] = lambdas
df["positions"] = positions
#df.to_csv("../Results/mean_dispersal_distances.csv")



df = pd.read_csv("../Results/mean_dispersal_distances.csv")

lambdas = df["lambdas"]
positions = df["positions"]

for y0 in [2,5,10,15,20]:
    # where pop > threshold
    x1_ = np.where(positions < y0)[0][0]
    x2_ = np.where(positions > y0)[0][-1]     
    y1 = positions[x1_]
    y2 = positions[x2_]
    x1 = lambdas[x1_]
    x2 = lambdas[x2_]

    # y = positions
    # x = lambdas

    # exponential 
    k = (np.log(y2) - np.log(y1)) / (x1 - x2)
    a = k * x1 + np.log(y1)
    x0 = (-np.log(y0) + a) / k


    Kernel = kernel_neg_exp_new(2**10, x0)
    #Kernel = cut_kernel(Kernel)
    dists = distances(2**10)
    position = np.sum(Kernel * dists)
    print(y0)
    print(x0)
    print(position)


#2 # desired distance
#1.1571953034774167 # lambda
#1.9999999996312976 # actual distance given by lambda

#5
#0.40940734064736645
#4.999999994609643

#10
#0.20120670217984254
#9.999999930553725

#15
#0.13369555710396816
#14.999999768315822

#20
#0.10015355023632337
#2

x0 = 1.158

Kernel = kernel_neg_exp_new(2**10, x0)
#Kernel = cut_kernel(Kernel)
dists = distances(2**10)
position = np.sum(Kernel * dists)

print(position)


