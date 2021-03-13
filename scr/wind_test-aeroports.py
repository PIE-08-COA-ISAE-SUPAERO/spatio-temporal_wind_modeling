#%% DATA Gathering
import TurboMaxWindElitePro as wind
import numpy as np
import matplotlib.pyplot as plt

folder = "G:/Mon Drive/PIE COA 08/Codes/TESTS/FICHIERS TESTS/AEROPORTS/"
aeroport = 'LFBO'

d = wind.wind()


d.create_wind_cube(folder+aeroport, aeroport)




#%% Test

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

