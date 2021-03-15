#%% Import
import TurboMaxWindElitePro as wind
import numpy as np
import matplotlib.pyplot as plt

#%% Data gathering*
folder = "C:/Users/loic-/OneDrive/Bureau/AEROPORTS/"
aeroport = ["LFRQ", "LFST", "LFRV"]

d = wind.wind()

for i in aeroport:
    d.create_wind_cube(folder + i, i)

#%% Test