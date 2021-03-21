#%% Import
import wind_class as wind
import numpy as np
import matplotlib.pyplot as plt
import os 
import csv
import xarray as xr
import seaborn as sns

import time 

# Data gathering
folder = "C:/Users/loic-/OneDrive/Bureau/TPLTE/"
# aeroport = ["LFBO", "LFKJ", "LFLP", "LFMY", "LFPB", "LFRC", "LFRQ", "LFRV", "LFST", "SOCA"]
aeroport = ["LFLP", "LFMY", "LFPB", "LFRC", "LFRQ", "LFRV", "LFST", "SOCA"]

debut = time.clock()
temps = []

wind_object = wind.wind()

# Wind creation
for i in aeroport:
    debut_int = time.clock()

    subfolder = folder  + i + '/'
    wind_object.create_wind_cube(subfolder, i)
    fin_int = time.clock()
    temps.append(fin_int - debut_int)
    print(i, ' : ', fin_int - debut_int)

fin = time.clock()
print(fin - debut, np.average(temps))

