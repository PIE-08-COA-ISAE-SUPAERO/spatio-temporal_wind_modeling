# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 10:48:01 2020

@author: Sébastien
"""

import numpy as np
import numpy.linalg as alg
import matplotlib.pyplot as plt
import xlrd
import numpy.random as rd
from cmath import *
from mpl_toolkits.mplot3d import axes3d 
import time


timedep = time.time()

#%% Génération d'un champ de vent fictif (décroissant selon z)
""" 
Le profil de vent U,V,W proposé ci-dessous sera utile pour faire les tests. 
Attention, on suppose ici que U, V, W désignent les trois composantes dans le 
même repère terrestre (et donc pas de notion de direction principale)
"""

X = np.linspace(0,400,10)
Y = np.linspace(0,400,10)
Z = np.linspace(0,400,10)
XX,YY,ZZ = np.meshgrid(X,Y,Z)

# U = 10+0*XX+0*YY+2/10*ZZ
# V = 0+0*XX+0*YY+0/10*ZZ
# W = 0+0*XX+0*YY+0/10*ZZ

U = 10+0*XX+0*YY+20*np.log(ZZ/10+1)
V = 0+0*XX+0*YY+0/10*ZZ
W = 0+0*XX+0*YY+0/10*ZZ

fig = plt.figure() 
ax = fig.gca(projection='3d') 
ax.quiver(XX, YY, ZZ, U, V, W, length=0.25, color = "blue") 
plt.show() 


#%% Génération d'un spectre localisé de turbulence


def Suu(f, vitesse, z) :
    
    lu = z/((0.177+0.00823*z)**(1.2))
    Umoy = (vitesse[0]**2+vitesse[1]**2+vitesse[2]**2)**(1/2)
    S = 4*lu/Umoy * 1/((1 + 70.7*(f*lu/Umoy)**2)**(5/6))
    
    return(S)


def Suv(f, vitesse, z) :
    
    lv = z/((0.177+0.00823*z)**(1.2))
    Umoy = (vitesse[0]**2+vitesse[1]**2+vitesse[2]**2)**(1/2)
    S = 4*lv/Umoy * (1 + 188.4*(2*f*lv/Umoy)**2)/((1 + 70.7*(2*f*lv/Umoy)**2)**(11/6))
    
    return(S)


def Suw(f, vitesse, z) :
    
    lw = z
    Umoy = (vitesse[0]**2+vitesse[1]**2+vitesse[2]**2)**(1/2)
    S = 4*lw/Umoy * (1 + 188.4*(2*f*lw/Umoy)**2)/((1 + 70.7*(2*f*lw/Umoy)**2)**(11/6))
    
    return(S)


vitesse_essai = [20,0,0]
hauteur_essai = 70
list_freq = np.linspace(0,10,100000)
list_Suu = [Suu(f, vitesse_essai, hauteur_essai) for f in list_freq]
list_Suv = [Suv(f, vitesse_essai, hauteur_essai) for f in list_freq]
list_Suw = [Suw(f, vitesse_essai, hauteur_essai) for f in list_freq]

plt.figure()
plt.plot(list_freq, list_Suu, label = "$S_u$")
plt.plot(list_freq, list_Suv, label = "$S_v$")
plt.plot(list_freq, list_Suw, label = "$S_w$")
plt.xscale("log")
plt.yscale("log")
plt.grid(True, which = "both", ls = "--", color= "0.75")
plt.xlabel("Fréquence $f$ (Hz)")
plt.ylabel("Densité spectrale de puissance")
plt.title("Spectre de turbulence")
plt.legend()


def freqTotime(vitesse, z, N, T, graph = 0): #retourne
    
    su, sv, sw = 1, 1, 1 #écarts types à modifier
    fs = N/T
    Xu, Xv, Xw = [], [], [] #variables aléatoires
    Tu, Tv, Tw = [], [], [] #turbulences
    temps = []
    tempsmax = int(T)
    
    for k in range(N):
        
        fk = k/N * fs
        suk = np.sqrt(T/(2*np.pi)*su**2*Suu(fk, vitesse, z))
        svk = np.sqrt(T/(2*np.pi)*sv**2*Suv(fk, vitesse, z))
        swk = np.sqrt(T/(2*np.pi)*sw**2*Suw(fk, vitesse, z))
        
        xu = rd.normal(0,suk)
        xv = rd.normal(0,svk)
        xw = rd.normal(0,swk)
        
        Xu.append(xu)
        Xv.append(xv)
        Xw.append(xw)
    
    for t in range(tempsmax) : #temps à définir intelligemment en fonction de l'échantillonage
        
        temps.append(t)
        Tut, Tvt, Twt = (vitesse[0]**2+vitesse[1]**2+vitesse[2]**2)**(1/2), 0, 0
        for k in range(N):
            Tut = Tut + 2*np.pi*fs/N * Xu[k] * np.cos(2*np.pi*k/N*t)
            Tvt = Tvt + 2*np.pi*fs/N * Xv[k] * np.cos(2*np.pi*k/N*t)
            Twt = Twt + 2*np.pi*fs/N * Xw[k] * np.cos(2*np.pi*k/N*t)
        Tu.append(Tut)
        Tv.append(Tvt)
        Tw.append(Twt)
    
    if graph == 1 :
        plt.figure()
        plt.plot(temps, Tu, label = "perturbation principale $\\bar{U}+\\tilde{u}$")
        plt.plot(temps, Tv, label = "perturbation latérale $\\tilde{v}$")
        plt.plot(temps, Tw, label = "perturbation verticale $\\tilde{w}$")
        plt.xlabel("Temps (s)")
        plt.ylabel("Vitesse du vent (m/s)")
        plt.title("Profil temporel du vent turbulent sur 10 minutes")
        plt.grid()
        plt.legend()
    
    return(Xu, Xv, Xw, Tu, Tv, Tw)


T_essai = 600
N_essai = 1000
Xu, Xv, Xw, Tu, Tv, Tw = freqTotime(vitesse_essai, hauteur_essai, N_essai, T_essai, 1)


#%% Perturbation sur le cube

timeint= time.time()

T_essai = 2
N_essai = 30

Ut = U.copy()
Vt = V.copy()
Wt = W.copy()


for x in range(np.shape(XX)[0]) :
    for y in range(np.shape(YY)[0]) :
        for z in range(np.shape(ZZ)[0]) :
            u, v, w = U[x,y,z], V[x,y,z], W[x,y,z]
            vitesse = [u, v, w]
            hauteur = ZZ[x,y,z]
            Xu, Xv, Xw, Tu, Tv, Tw = freqTotime(vitesse, hauteur, N_essai, T_essai)
            Ut[x,y,z] = Tu[1]
            Vt[x,y,z] = Tv[1]
            Wt[x,y,z] = Tw[1]


print("tailles maillage/stationnaire/turbulent", np.shape(XX),np.shape(U), np.shape(Ut))


fig = plt.figure() 
ax = fig.gca(projection='3d') 
ax.quiver(XX, YY, ZZ, U, V, W, length = 0.35, color = "limegreen", label = "Stationnaire") 
ax.quiver(XX, YY, ZZ, Ut, Vt, Wt, length = 0.35, color = "red", label = "Turbulent")
ax.legend() 
plt.show() 


print(U[1,2,3], "versus", Ut[1,2,3])
print(V[1,2,3], "versus", Vt[1,2,3])
print(W[1,2,3], "versus", Wt[1,2,3])

timefin = time.time()

print("total", timefin-timedep, "partiel", timefin-timeint)



"""
Une question qui se pose : après avoir inséré de la turbulence en tout point
du maillage, on ne respecte a priori plus les équations de la dynamique (chaque
perturbation, aléatoire, ne dépend que de la hauteur au mieux). Faut-il garder
alors le vent après la perturbation ou doit-on envisager de remailler (en faisant 
par exemple passer une conservation de la masse ou du lissage) ?

10 -> 15 s  (3)
20 -> 36 s  (21)  (x2, x8)
40 -> 131 s (118) (x2, x5)
80 -> 864 s (838) (x2, x7)
Opération sur du O(^3)

Doit être considéré pour être plus précis:
    - la détermination du nombre de composantes N à prendre
    - l'assurance que la direction est bonne (changement de base)
    - un essai sur un vent stationnaire effectif

"""


















#%% Inspiration grilles 3D

# fig = plt.figure() 
# ax = fig.gca(projection='3d') 

# x, y, z = np.meshgrid(np.arange(-0.8, 1, 0.2), 
#          np.arange(-0.8, 1, 0.2), 
#          np.arange(-0.8, 1, 0.8)) 

# u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z) 
# v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z) 
# w = (np.sqrt(2.0/3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) * 
#     np.sin(np.pi * z)) 

# print(type(u))
# print(type(x))
# print(np.shape(x),np.shape(u))

# ax.quiver(x, y, z, u, v, w, length=0.1) 

# plt.show() 