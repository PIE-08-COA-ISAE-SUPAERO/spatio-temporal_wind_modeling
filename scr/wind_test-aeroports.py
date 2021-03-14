#%% Import
import TurboMaxWindElitePro as wind
import numpy as np
import matplotlib.pyplot as plt

#%% Data gathering*
folder = "G:/Mon Drive/PIE COA 08/Codes/TESTS/FICHIERS TESTS/AEROPORTS/"
aeroport = 'LFBO'

d = wind.wind()

d.create_wind_cube(folder+aeroport, aeroport)

#%% Test