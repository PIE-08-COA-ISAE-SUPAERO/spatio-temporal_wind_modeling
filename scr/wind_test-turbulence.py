#%% DATA Gathering
import wind_class as wind
import numpy as np
import matplotlib.pyplot as plt

import xlrd 
import seaborn as sns
import os

#%%
ROOT = 'G:/Mon Drive/PIE COA 08/Codes/PGM/spatio-temporal_wind_modeling/spatio-temporal_wind_modeling/data/Test_phase1_doctorante/'
DT = 0.1

def fichier_mat(altitude, short = True):
    path = ROOT + 'INPUT DATA/mat_{}m'.format(altitude)
    if short : path += '_short'
    return path + '.xlsx'

def fichier_data_wn(hours):
    return ROOT + 'OUTPUT DATA/H{}/exported_data_.json'.format(hours)

coord_station = (43.1498753, 0.3560498)
liste_altitudes = [30, 45, 60]

d = wind.wind()

def get_data(hours, altitude, duration_s, short = True) :    
    print('Begin')
    # Field data
    input_data_wb = xlrd.open_workbook(fichier_mat(altitude, short))
    input_data_ws = input_data_wb.sheet_by_index(0)
    print('File imported')
    heure_ws = hours * 3600 * 1000
    n = 10 * duration_s

    u_ws = np.zeros(n)
    v_ws = np.zeros(n)
    w_ws = np.zeros(n)

    #We look for the correct row and collect the value
    for i in range(input_data_ws.nrows):
        if input_data_ws.cell_value(i, 1) == heure_ws :
            print('Hour found')
            for j in range(n) :
                u_ws[j] = input_data_ws.cell_value(i+j, 2)
                v_ws[j] = input_data_ws.cell_value(i+j, 3)
                w_ws[j] = input_data_ws.cell_value(i+j, 4)

    #We get the simulated data
    print('Get simulation')
    d.import_wind_cube(fichier_data_wn(hours))

    print('Get turbulence')
    _, _, _, u_wn, v_wn, w_wn = d.profil_turbulence(coord_station[0], coord_station[1], altitude, DT, plot=False)

    print('End')
    return np.array([u_wn, u_ws]), np.array([v_wn, v_ws]), np.array([w_wn, w_ws])

H = 1
alt = liste_altitudes[1]
duration = 15*60
u, v, w = get_data(H, alt, duration, short= False)
#%%

def print_courbe(data, name, plot = True):
    simule, reelle = data

    n = len(simule)
    x = np.linspace(0, 15, n)    

    moyenne_sim = [np.average(simule)] * n
    moyenne_relle = [np.average(reelle)] * n

    plt.figure()
    plt.plot(x, moyenne_sim, 'r')
    plt.plot(x, simule, '--r', label = 'simulé', )

    plt.plot(x, moyenne_relle, 'b')
    plt.plot(x, reelle, '--b', label = 'reele')
    plt.legend()
    plt.show()

def print_comparaison(data, name, plot = True):
    simule, reelle = data

    plt.figure()
    plt.plot(reelle, simule)
    plt.xlabel('reel')
    plt.ylabel('simule')
    plt.show()

print_courbe(u, 'u')
print_courbe(v, 'v')
print_courbe(w, 'w')