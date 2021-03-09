import TurboMaxWindElitePro as wind
import numpy as np

folder = 'G:/Mon Drive/PIE COA 08/Codes/TESTS/FICHIERS TESTS/Données doctorante/OUPUT DATA/H'
coord_station = (43.1498753, 0.3560498)

N = 1

d = wind.wind()

# for i in range(N):
#    folder_i = folder + str(i+1)
#    d.create_wind_cube(folder_i, folder_i)

d.import_wind_cube(folder + str(N) + '/exported_data_.json')

# Wind_speed = d._wind_cube["Wind_speed"]

# U = Wind_speed[:,0]
# V = Wind_speed[:,1]
# W = Wind_speed[:,2]

# u_moy =np.average(U)
# v_moy =np.average(V)
# w_moy =np.average(W)

# print(u_moy, v_moy, w_moy)
# print(d.get_point(coord_station[0], coord_station[1], elevation= 30, plot = True))
d.plot_wind_surface("x", coord_station, 30, plot = True)