# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 14:56:01 2020

@author: Q.Abeille (2020)
"""
#%% Imports
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('presentation.mplstyle')
from scipy import interpolate as intp
from scipy import optimize as opt

#%% Recupération des données

points = np.loadtxt('test_wind.vtk', skiprows=7, max_rows=78840)                # Positions des points de calcul
vectors = np.loadtxt('test_wind.vtk', skiprows=78850)                           # Champ de vitesse en chaque point

#%% Traitement des données
nb_points = np.size(points,0)
# Décomposition selon les 3 axes cartésiens
X = np.array(points[:,0])
Y = np.array(points[:,1])
Z = np.array(points[:,2])

# Récupération des composantes du vent
U = np.array(vectors[:,0])
V = np.array(vectors[:,1])
W = np.array(vectors[:,2])

# Maillage
X_tick = np.unique(X)
Y_tick = np.unique(Y)

# Champ de vitesse pour chaque composante du vent
U_field = np.array([X,Y,Z,U])
U_field = np.transpose(U_field)

V_field = np.array([X,Y,Z,V])
V_field = np.transpose(V_field)

W_field = np.array([X,Y,Z,W])
W_field = np.transpose(W_field)

#%% Récupération du vent en fonction de l'altitude en un point donné
x_pos = X_tick[0]
y_pos = Y_tick[0]

def get_vert_profile(x_pos,y_pos):
    
    cond_X = U_field[:,0]==x_pos                                                # Elts situés en x_pos
    cond_Y = U_field[:,1]==y_pos                                                # Elts situés en y_pos
    cond = np.multiply(cond_X, cond_Y)                                          # Elts situés en x_pos, y_pos
    
    ZUarg_pos = np.argwhere(cond==True)                                         # Récupération de la position des elts
    
    # Tableaux de sortie
    Zxy = np.array([])
    Uxy = np.array([])
    Vxy = np.array([])
    Wxy = np.array([])
    
    # Remplissage des tableaux
    for k in range (0,len(ZUarg_pos)):
        cur_arg = ZUarg_pos[k]
        Zxy = np.append(Zxy,U_field[cur_arg,2])
        Uxy = np.append(Uxy,U_field[cur_arg,3])
        Vxy = np.append(Vxy,V_field[cur_arg,3])
        Wxy = np.append(Wxy,W_field[cur_arg,3])
    return(Zxy,Uxy,Vxy,Wxy)

[Zxy,Uxy,Vxy,Wxy] = get_vert_profile(x_pos,y_pos)

# plt.figure()
# plt.plot(Uxy,Zxy,'b-o',label='U')
# plt.grid()

# plt.figure()
# plt.plot(Vxy,Zxy,'g-o',label='V')
# plt.grid()

# plt.figure()
# plt.plot(Wxy,Zxy,'y-o',label='W')
# plt.grid()


#%% Interpolation sur U

# Loi puissance
def power_law(z,U10,alpha):
    return U10*np.power((z-1554.23)/10., alpha)

# Loi log
def log_law(z,U10,z0):
    return U10/np.log(10./0.05) * np.log((z-1554.23)/0.05)

start = 14          # indice de début pour la séléction des données
alt_max = 5564      # altitude max d'interpolation

# Création des tableaux pour l'interpolation puis l'extrapolation
Z_extrap=np.concatenate((Zxy[start:20],[alt_max]))
U_extrap = np.concatenate((Uxy[start:20],[0.]))
Z_stud = np.arange(Zxy[start],alt_max,50)

# Loi puissance
best_vals1, covar1 = opt.curve_fit(power_law, Z_extrap, U_extrap, p0=[2.5,0.128])
u_fit1 = power_law(Z_stud,best_vals1[0],best_vals1[1])

# Loi log
best_vals2, covar2 = opt.curve_fit(log_law, Z_extrap, U_extrap, p0=[2.5,0.05])
u_fit2 = log_law(Z_stud,best_vals2[0],best_vals2[1])

# Loi polynomiale degré 7
u_params = np.polyfit(Z_extrap, U_extrap, 7)
u_fit3 = np.poly1d(u_params)

# Interpolation linéaire
f_linearU = intp.interp1d(Z_extrap,U_extrap,kind='slinear')

# Interpolation cubique
f_cubicU = intp.interpolate.interp1d(Z_extrap,U_extrap,kind='cubic')

# Présentation des résultats
plt.figure('U',figsize=(14,10))
plt.plot(Uxy,Zxy,'b-o',label='U')
plt.plot(u_fit1,Z_stud,'r:x',label='power_law')
plt.plot(u_fit2,Z_stud,'g:x',label='log_law')
#plt.plot(u_fit3(Z_stud),Z_stud,'y:x',label='poly7')
plt.plot(f_linearU(Z_stud),Z_stud,'m:x',label = 'interp1d_slinear')
#plt.plot(f_cubicU(Z_stud),Z_stud,'b:x',label = 'interp1d_cubic')
plt.legend(fontsize=16)
plt.xlabel('Vitesse algébrique (km/h)')
plt.ylabel('altitude (m)')
plt.grid()

#%% Interpolation sur V

# Loi puissance
def power_law(z,V10,alpha):
    return V10*np.power((z-1554.23)/10., alpha)

# Loi log
def log_law(z,V10,z0):
    return V10/np.log(10./0.05) * np.log((z-1554.23)/0.05)

# Création des tableaux pour interpolation puis extrapolation
Z_extrap=Zxy[start:20]
V_extrap = Vxy[start:20]
Z_stud = np.arange(Zxy[start],alt_max,50)

# Loi puissance
best_vals1, covar1 = opt.curve_fit(power_law, Z_extrap, V_extrap, p0=[2.5,0.128])
v_fit1 = power_law(Z_stud,best_vals1[0],best_vals1[1])

# Loi log
best_vals2, covar2 = opt.curve_fit(log_law, Z_extrap, V_extrap, p0=[2.5,0.05])
v_fit2 = log_law(Z_stud,best_vals2[0],best_vals2[1])

# Loi poly degré 7
v_params = np.polyfit(Z_extrap, V_extrap, 7)
v_fit3 = np.poly1d(v_params)

# Interpolation linéaire
f_linearV = intp.interp1d(Z_extrap,V_extrap,kind='slinear',bounds_error=False,fill_value='extrapolate')

# Interpolation cubique
f_cubicV = intp.interpolate.interp1d(Z_extrap,V_extrap,kind='cubic',bounds_error=False, fill_value = 'extrapolate')

# Présentation des résultats
plt.figure('V',figsize=(14,10))
plt.plot(Vxy,Zxy,'b-o',label='V')
plt.plot(v_fit1,Z_stud,'r:x',label='power_law')
plt.plot(v_fit2,Z_stud,'g:x',label='log_law')
#plt.plot(v_fit3(Z_stud),Z_stud,'y:x',label='poly7')
plt.plot(f_linearV(Z_stud),Z_stud,'m:x',label = 'interp1d_slinear')
#plt.plot(f_cubicV(Z_stud),Z_stud,'b:x',label = 'interp1d_cubic')
plt.legend(fontsize=16)
plt.xlabel('Vitesse algébrique (km/h)')
plt.ylabel('altitude (m)')
plt.grid()

#%% Interpolation sur W

# Loi puissance
def power_law(z,W10,alpha):
    return W10*np.power((z-1554.23)/10., alpha)

# Loi log
def log_law(z,W10,z0):
    return W10/np.log(10./0.05) * np.log((z-1554.23)/0.05)

# Création des tableaux pour interpolation puis extrapolation
Z_extrap=np.concatenate((Zxy[start:20],[alt_max]))
W_extrap = np.concatenate((Wxy[start:20],[0.]))
Z_stud = np.arange(Zxy[start],alt_max,50)

# Loi puissance
best_vals1, covar1 = opt.curve_fit(power_law, Z_extrap, W_extrap, p0=[2.5,0.128])
w_fit1 = power_law(Z_stud,best_vals1[0],best_vals1[1])

# Loi log
best_vals2, covar2 = opt.curve_fit(log_law, Z_extrap, W_extrap, p0=[2.5,0.05])
w_fit2 = log_law(Z_stud,best_vals2[0],best_vals2[1])

# Loi poly degré 7
w_params = np.polyfit(Z_extrap, W_extrap, 7)
w_fit3 = np.poly1d(w_params)

# Interpolation linéaire
f_linearW = intp.interp1d(Z_extrap,W_extrap,kind='slinear',bounds_error=False, fill_value = 'extrapolate')

# Interpolation cubique
f_cubicW = intp.interpolate.interp1d(Z_extrap,W_extrap,kind='cubic')

# Présentation des résultats
plt.figure('W',figsize=(14,10))
plt.plot(Wxy,Zxy,'b-o',label='W')
plt.plot(w_fit1,Z_stud,'r:x',label='power_law')
plt.plot(w_fit2,Z_stud,'g:x',label='log_law')
#plt.plot(w_fit3(Z_stud),Z_stud,'y:x',label='poly7')
plt.plot(f_linearW(Z_stud),Z_stud,'m:x',label = 'interp1d_slinear')
#plt.plot(f_cubicW(Z_stud),Z_stud,'b:x',label = 'interp1d_cubic')
plt.legend(fontsize=16)
plt.xlabel('Vitesse algébrique (km/h)')
plt.ylabel('altitude (m)')
plt.grid()

#%% Norme du vent

wind_norm =[]

# Calcul de la norme
for k in range(len(Uxy)):
    norm_z = np.sqrt(Uxy[k]**2 + Vxy[k]**2 + Wxy[k]**2)
    wind_norm = np.append(wind_norm,norm_z)
    
# Loi puissance
def power_law(z,wind10,alpha):
    return wind10*np.power((z-1554.23)/10., alpha)

# Loi :og
def log_law(z,wind10,z0):
    return wind10/np.log(10./0.05) * np.log((z-1554.23)/0.05)

# Tableaux pour interpolation puis extrapolation
Z_extrap= Zxy[start:20]
wind_extrap = wind_norm[start:20]

# Loi puissance
best_vals1, covar1 = opt.curve_fit(power_law, Z_extrap, wind_extrap, p0=[2.5,0.128])
wind_fit1 = power_law(Zxy,best_vals1[0],best_vals1[1])

# Loi log
best_vals2, covar2 = opt.curve_fit(log_law, Z_extrap, wind_extrap, p0=[2.5,0.05])
wind_fit2 = log_law(Zxy,best_vals2[0],best_vals2[1])

# Interpolation linéaire
f_linearNorm = intp.interp1d(Z_extrap,wind_extrap,kind='slinear',bounds_error=False, fill_value = 'extrapolate')

# Présentation des résultats
plt.figure('Norme',figsize=(14,10))
plt.plot(wind_norm,Zxy,'b-o', label='wind_norm')
plt.plot(wind_fit1,Zxy,'r:x', label='power_law')
plt.plot(wind_fit2,Zxy,'g:x', label='log_law')
plt.plot(f_linearNorm(Zxy),Zxy, 'm:x', label='interp1d_slinear')
plt.grid()
plt.xlabel('Vitesse (km/h)')
plt.ylabel('altitude (m)')
plt.legend(fontsize=16)