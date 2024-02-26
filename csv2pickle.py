#!/usr/bin/env python

"""
take csv and save each row as a pickle
"""

import pandas as pd 
import pickle 

df = pd.read_csv("../Data/inv_speed_6.536_G_sim_params.csv")

for i in range(df.shape[0]):
    #i=0
    # subset row and convert to dict 
    dic = df.iloc[[i]].to_dict("records")[0]
    # save dic as pickle
    outfilename = "../Data/input" + str(i) + ".pkl"
    with open(outfilename, 'wb') as f:
        pickle.dump(dic, f)

#with open("test_parallel/output0.pkl", 'rb') as f:
#    A = pickle.load(f)
