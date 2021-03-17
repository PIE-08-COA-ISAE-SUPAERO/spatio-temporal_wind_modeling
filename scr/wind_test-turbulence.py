#%% DATA Gathering
import wind_class as wind
import numpy as np
import matplotlib.pyplot as plt

import xlrd 
import seaborn as sns
import os

#%%
ROOT = 'G:/Mon Drive/PIE COA 08/Codes/PGM/spatio-temporal_wind_modeling/spatio-temporal_wind_modeling/data/Test_phase1_doctorante/'

def fichier_mat(altitude, short = True):
    path = ROOT + 'INPUT DATA/mat_{}m'.format(altitude)
    if short : path += '_short'
    return path + '.xlsx'

def fichier_data_wn(hours):
    return ROOT + 'OUTPUT DATA/H{}/exported_data_.json'.format(hours)

coord_station = (43.1498753, 0.3560498)
liste_altitudes = [30, 45, 60]

d = wind.wind()

h = 1
i = 1
altitude = liste_altitudes[i]

input_data_wb = xlrd.open_workbook(fichier_mat(altitude, short = True))
input_data_ws = input_data_wb.sheet_by_index(0)

#On recupere la ligne de l'heure
#On fait 5min a 10ms de pas et on recupèrere toutes les donneés, on les compare

#Au pire, on lance la simulation N fois et on moyenne