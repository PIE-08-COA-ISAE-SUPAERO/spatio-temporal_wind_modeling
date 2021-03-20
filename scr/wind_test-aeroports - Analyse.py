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
ROOT = "G:/Mon Drive/PIE COA 08/Codes/PGM/spatio-temporal_wind_modeling/spatio-temporal_wind_modeling/data/test_phase2_aeroports/"

liste_aeroport = ["LFBO", "LFKJ", "LFLP", "LFMY", "LFPB", "LFRC", "LFRQ", "LFRV", "LFST", "SOCA"]
liste_folder = ["2021-03-17_1500", "2021-03-19_1800"]

fichier_aeroport = 'Liste_aeroports.csv'

N_folder = len(liste_folder)
i_aeroport = 0
N_aeroport = len(liste_aeroport)

wind_speed = []
wind_direction = []

wind_object = wind.wind()

for i in range(N_folder):
    folder = ROOT + liste_folder[i] + '/'
    i_aeroport = 0

    with open(folder + fichier_aeroport) as csvfile :
        spamreader = csv.reader(csvfile, delimiter = ';')

        for row in spamreader:
            if row[1] != 'Nom' :
                stid = row[0]

                subfolder = folder + stid + '/'

                #Get the airport data
                lon_aeroport = float(row[2])
                lat_aeroport = float(row[3])
                alt_aeroport = float(row[4])
                wind_speed_aeroport = float(row[6])
                wind_direction_aeroport = float(row[7]) if row[7] != '' else 0

                #Get the wn data
                files = wind.file_list_by_extension(subfolder, '.json')
                data_file = [_ for _ in files if 'exported' in _]
                
                wind_object.import_wind_cube(subfolder + data_file[0])
                wn_point_param = wind_object.get_point(lat_aeroport, lon_aeroport, elevation = alt_aeroport, plot = False)
                
                wind_speed.append([wn_point_param[-3], wind_speed_aeroport])
                if row[7] != '' : wind_direction.append([wn_point_param[-1], wind_direction_aeroport])
                
                print(row, (i_aeroport+1) / N_aeroport * 100, '%')
                i_aeroport += 1

    print('---------------', liste_folder[i], (i+1) / N_folder* 100, '%')

print('########','Done', sep = '\n')

#%% Tests
def hist_err(data_brute, name = ''):
    data_ecart = [(x[0] - x[1]) for x in data_brute]
    N = int(len(data_ecart)*0.6)
    
    sns.distplot(data_ecart, bins = N)
    plt.ylabel('Probalbilité %')
    plt.xlabel('Erreur commise (m/s)')
    plt.show()

def err_droite(data, name, plot = True):
    data = np.transpose(np.array(data)) #[[sim], [reel]]
    sim, reel = data

    # a, _, _, _ =  np.linalg.lstsq(reel[:,np.newaxis], sim)
    r = np.corrcoef(data)

    a = 1
    x = np.linspace(min(0,np.min(reel)), np.max(reel))
    y = a * x

    print(name,' : ', r[0,1])
    if plot :
        plt.figure(name)
        
        plt.plot(sim, reel, '.')
        plt.plot(x,y)
        
        plt.xlabel('Réel (m/s)')
        plt.ylabel('Simulé (m/s)')

        plt.title('Variable : ${}$, $a$ = {} (forcé), $r^2$ = {}'.format(name, int(a*1000)/1000, int(r[0,1]**2*1000)/1000))
        
        plt.legend()
        plt.show()

def print_err_simu(data, name, plot = True):
    
    data = np.transpose(np.array(data)) #[[sim], [reel]]
    sim, reel = data

    a, _, _, _ =  np.linalg.lstsq(reel[:,np.newaxis], sim)
    r = np.corrcoef(data)

    a = a[0]
    x = np.linspace(min(0,np.min(reel)), np.max(reel))
    y = a * x

    print(name,' : ', r[0,1])
    if plot :
        plt.figure(name)
        
        plt.plot(sim, reel, '.')
        plt.plot(x,y)
        
        plt.xlabel('Réel (m/s)')
        plt.ylabel('Simulé (m/s)')

        plt.title('Variable : ${}$, $a$ = {}, $r^2$ = {}'.format(name, int(a*1000)/1000, int(r[0,1]**2*1000)/1000))
        
        plt.legend()
        plt.show()

err_droite(wind_speed, 'Vitesse', True)
err_droite(wind_direction, 'Direction', True)
print_err_simu(wind_speed, 'Vitesse', True)
print_err_simu(wind_direction, 'Direction', True)
hist_err(wind_speed, 'vitesse')