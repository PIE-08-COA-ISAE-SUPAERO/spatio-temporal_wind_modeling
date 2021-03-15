#%% Import
import TurboMaxWindElitePro as wind
import numpy as np
import matplotlib.pyplot as plt

#%% Data gathering*
folder = "C:/Users/abeil/source/repos/spatio-temporal_wind_modeling/data/Tests_airports/"
aeroport = "LFBO"

d = wind.wind()

d.create_wind_cube(folder+aeroport, aeroport, "")

#%% Test