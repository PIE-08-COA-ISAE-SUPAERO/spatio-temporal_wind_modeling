#%% DATA Gathering
import TurboMaxWindElitePro as wind
import numpy as np
import matplotlib.pyplot as plt

# from scipy.stats import linregress

import xlrd 
import seaborn as sns

#%%
folder = 'G:/Mon Drive/PIE COA 08/Codes/TESTS/FICHIERS TESTS/Données doctorante/OUPUT DATA/H'
coord_station = (43.1498753, 0.3560498)

input_data_wb = xlrd.open_workbook('G:/Mon Drive/PIE COA 08/Codes/TESTS/FICHIERS TESTS/Données doctorante/INPUT DATA/data_test_doctorante.xls')

altitude_station = [30, 45, 60]
N = 20

u = np.zeros((N, 3, 2)) #h, alt, (0 = simule, 1 = reel)
v = np.zeros((N, 3, 2))
w = np.zeros((N, 3, 2))
norme = np.zeros((N, 3, 2))
norme_plan = np.zeros((N, 3, 2))

d = wind.wind()

for i in range(N):
    folder_i = folder + str(i+1)
    d.import_wind_cube(folder_i + '/exported_data_.json')

    for j in range(3):
            input_data_ws = input_data_wb.sheet_by_index(j)

            u_,v_,w_,norme_, norme_plan_,_ = d.get_point(coord_station[0], coord_station[1], elevation= altitude_station[j], plot = False)
            u[i, j, 0] = u_
            v[i, j, 0] = v_
            w[i, j, 0] = w_
            norme[i, j, 0] = norme_
            norme_plan[i, j, 0] = norme_plan_

            u_mes = input_data_ws.cell_value(i+1,1)
            v_mes = input_data_ws.cell_value(i+1,2)
            w_mes = input_data_ws.cell_value(i+1,3)
            norme_mes = input_data_ws.cell_value(i+1,4)
            norme_plan_mes = np.sqrt(u_mes**2 + v_mes**2)

            u[i, j, 1] = u_mes
            v[i, j, 1] = v_mes
            w[i, j, 1] = w_mes
            norme[i, j, 1] = norme_mes
            norme_plan[i, j, 1] = norme_plan_mes
 
    print((i+1)/N*100, '%')

#%% Reshaping of the data
U_Fin = np.zeros((3*N, 2))
V_Fin = np.zeros((3*N, 2))
W_Fin = np.zeros((3*N, 2))
Norme_Fin = np.zeros((3*N, 2))
Norme_Plan_Fin = np.zeros((3*N, 2))

i_heure = i_alt = i_glob = 0

while i_heure < N :
    while i_alt < 3 :
        U_Fin[i_glob,:] = u[i_heure, i_alt, :]
        V_Fin[i_glob,:] = v[i_heure, i_alt, :]
        W_Fin[i_glob,:] = w[i_heure, i_alt, :]
        Norme_Fin[i_glob,:] = norme[i_heure, i_alt, :]
        Norme_Plan_Fin[i_glob,:] = norme_plan[i_heure, i_alt, :]

        i_alt += 1
        i_glob += 1
    i_alt =0
    i_heure += 1

U_err = U_Fin[:,1] - U_Fin[:,0]
V_err = V_Fin[:,1] - V_Fin[:,0]
W_err = W_Fin[:,1] - W_Fin[:,0]
Norme_err = Norme_Fin[:,1] - Norme_Fin[:,0]
Norme_Plan_Err = Norme_Plan_Fin[:,1] - Norme_Plan_Fin[:,0]

#%% Test 

def hist_err(data, name = ''):
    data = np.ravel(data)
    # N = int(len(data)*1)
    
    # print(N, len(data))

    # min_x = np.min(data)
    # max_x = np.max(data)

    # x = np.linspace(min_x, max_x, N)
    # y = np.zeros(N)

    # for i in range(N):
    #     for j in range(N-1):
    #         if x[j] <= data[i] and data[i] <= x[j+1]:
    #             y[j] += 1
    # y = y / len(data) *100

    # x_filtre = []
    # y_filtre = []

    # for i in range(N):
    #     if y[i] != 0:
    #         x_filtre.append(x[i])
    #         y_filtre.append(y[i])

    # plt.figure(name)
    # plt.plot(x_filtre, y_filtre, '-')

    # plt.xlabel('Erreur résiduelle (m/s)')
    # plt.ylabel("Probabilité d'erreur (%)")

    # plt.xlim([np.min(x_filtre), np.max(x_filtre)])
    # plt.ylim([0, np.max(y_filtre)*1.1])

    # plt.show()
    sns.distplot(data)
    
hist_err([U_err, V_err, W_err, Norme_err, Norme_Plan_Err])


def print_err_simu(data, name, plot = True):
        
    sim, reel = data[:,0],data[:,1]

    # if name == 'w' :
    #     sim = [i *1000 for i in sim]
    #     reel = [i*1000 for i in reel]

    a, _, _, _ =  np.linalg.lstsq(reel[:,np.newaxis], sim)

    a = a[0]
    x = np.linspace(min(0,np.min(reel)), np.max(reel))
    y = a * x

    r = np.corrcoef(np.transpose(data))

    if plot : 
        plt.figure(name)
        plt.plot(reel, sim, 'o')
        plt.plot(x,y)

        plt.xlim([min(0,np.min(reel)), np.max(reel)])
        plt.ylim([min(0,np.min(sim)),np.max(sim)])

        plt.xlabel('Réel (m/s)')
        plt.ylabel('Simulé (m/s)')

        plt.title('Variable : ${}$, $a$ = {}, $r^2$ = {}'.format(name, int(a*1000)/1000, int(r[0,1]**2*1000)/1000))
        
        plt.legend()
        plt.show()

# print_err_simu(U_Fin, 'u', True)
# print_err_simu(V_Fin, 'v', True)
# print_err_simu(W_Fin, 'w', True)
# print_err_simu(Norme_Fin, 'norme', True)
# print_err_simu(Norme_Plan_Fin, 'norme plane', True)

