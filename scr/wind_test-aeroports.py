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
folder = "G:/Mon Drive/PIE COA 08/Codes/PGM/spatio-temporal_wind_modeling/spatio-temporal_wind_modeling/data/test_phase2_aeroports/"
# folder = "C:/Users/loic-/OneDrive/Bureau/test_phase2_aeroports/"
liste_aeroport = folder  + 'Liste_aeroports.csv'

# debut = time.clock()
# temps = []

wind_object = wind.wind()

###################
# Wind creation
# for i in aeroport:
#     debut_int = time.clock()

#     subfolder = folder  + i + '/'
#     wind_object.create_wind_cube(subfolder, i)
#     fin_int = time.clock()
#     temps.append(fin_int - debut_int)
#     print(i, ' : ', fin_int - debut_int)

# fin = time.clock()
# print(fin - debut, np.average(temps))

###################
#Test

wind_speed = []
wind_planed_speed = []
wind_direction = []

with open(liste_aeroport) as csvfile :
    spamreader = csv.reader(csvfile, delimiter = ';')

    for row in spamreader:
        if row[1] != 'Nom':
            stid = row[0]

            subfolder = folder + stid + '/'

            print(row)

            lon_aeroport = float(row[2])
            lat_aeroport = float(row[3])
            alt_aeroport = float(row[4])
            wind_speed_aeroport = float(row[6])
            wind_direction_aeroport = float(row[7]) if row[7] != '' else 0

            files = wind.file_list_by_extension(subfolder, '.json')
            data_file = [_ for _ in files if 'exported' in _]

            wind_object.import_wind_cube(subfolder + data_file[0])
            wn_point_param = wind_object.get_point(lat_aeroport, lon_aeroport, elevation = alt_aeroport, plot = False)
            
            print(wn_point_param)
            wind_speed.append([wn_point_param[-3], wind_speed_aeroport])
            wind_planed_speed.append([wn_point_param[-2], wind_speed_aeroport])
            if row[7] != '' : wind_direction.append([wn_point_param[-1], wind_direction_aeroport])

print('Done')
#%% Tests
def hist_err(data_brute, name = ''):
    data_ecart = [(x[0] - x[1])/x[0]*100 for x in data_brute]
    N = int(len(data_ecart)*0.6)
    
    sns.distplot(data_ecart, bins = N)
    plt.ylabel('Probalbilité %')
    plt.xlabel('Erreur commise (%)')
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

# err_droite(wind_speed, 'Vitesse', True)
# err_droite(wind_planed_speed, 'Vitesse plane', True)
# print_err_simu(wind_speed, 'Vitesse', True)
# print_err_simu(wind_planed_speed, 'Vitesse plane', True)
# print_err_simu(wind_direction, 'Direction', True)
hist_err(wind_speed, 'vitesse')