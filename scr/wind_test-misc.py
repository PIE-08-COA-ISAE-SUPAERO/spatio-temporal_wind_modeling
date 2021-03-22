#%% Import
import wind_class as wind
import numpy as np
import matplotlib.pyplot as plt


# Data gathering
ROOT = "G:/Mon Drive/PIE COA 08/Codes/PGM/spatio-temporal_wind_modeling/spatio-temporal_wind_modeling/data/test_phase2_aeroports/"

wind_object = wind.wind()

wind_object.import_wind_cube(ROOT + 'exported_data_Missoula_2021-03-21.json')

x = wind_object.cube_coordinates()
lat = x["Latitude"]["max"]
lon = x["Longitude"]["max"]
print(wind.flat_distance_point(wind_object._location, [lat, lon]))
