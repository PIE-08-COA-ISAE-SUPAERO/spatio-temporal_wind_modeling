# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 12:03:27 2021

@author: Sébastien
"""


import numpy as np
import numpy.linalg as alg
from numpy.fft import fft
import matplotlib.pyplot as plt
import xlrd
from cmath import *
import time
from scipy.optimize import curve_fit, fsolve

essai = 0
Ndep = 50000

#%% LECTURE DES FICHIERS

if essai == 1 :
    fichier30 = 'mat_30m.xlsx'
    fichier45 = 'mat_45m.xlsx'
    fichier60 = 'mat_60m.xlsx'
else :
    fichier30 = 'mat_30m_short.xlsx'
    fichier45 = 'mat_45m_short.xlsx'
    fichier60 = 'mat_60m_short.xlsx'

lecture30 = xlrd.open_workbook(fichier30)
feuille30 = lecture30.sheet_by_index(0)
nlig = feuille30.nrows
ncol = feuille30.ncols
temps30, u30, v30, w30, norme30, direction30 = [], [], [], [], [], []
for k in range(1, nlig):
    temps_lu = feuille30.cell_value(k, 1)
    u_lu = feuille30.cell_value(k, 2)
    v_lu = feuille30.cell_value(k, 3)
    w_lu = feuille30.cell_value(k, 4)
    if u_lu != "NaN" and v_lu != "NaN" and w_lu != "NaN" :
        temps30.append(float(temps_lu)/1000) #temps en s
        u30.append(float(u_lu))
        v30.append(float(v_lu))
        w30.append(float(w_lu))
        norme30.append(np.sqrt(u30[-1]**2 + v30[-1]**2 + w30[-1]**2 ))
        if u30[-1]**2 > 0.001 : 
            direction30.append(np.arctan(v30[-1]/u30[-1])*180/np.pi)
        else :
            direction30.append(90)
        

lecture45 = xlrd.open_workbook(fichier45)
feuille45 = lecture45.sheet_by_index(0)
nlig = feuille45.nrows
ncol = feuille45.ncols
temps45, u45, v45, w45, norme45, direction45 = [], [], [], [], [], []
for k in range(1, nlig):
    temps_lu = feuille45.cell_value(k, 1)
    u_lu = feuille45.cell_value(k, 2)
    v_lu = feuille45.cell_value(k, 3)
    w_lu = feuille45.cell_value(k, 4)
    if u_lu != "NaN" and v_lu != "NaN" and w_lu != "NaN" :
        temps45.append(float(temps_lu)/1000) #temps en s
        u45.append(float(u_lu))
        v45.append(float(v_lu))
        w45.append(float(w_lu))
        norme45.append(np.sqrt(u45[-1]**2 + v45[-1]**2 + w45[-1]**2 ))
        if u45[-1]**2 > 0.001 : 
            direction45.append(np.arctan(v45[-1]/u45[-1])*180/np.pi)
        else :
            direction45.append(90)

lecture60 = xlrd.open_workbook(fichier60)
feuille60 = lecture60.sheet_by_index(0)
nlig = feuille60.nrows
ncol = feuille60.ncols
temps60, u60, v60, w60, norme60, direction60 = [], [], [], [], [], []
for k in range(1, nlig):
    temps_lu = feuille60.cell_value(k, 1)
    u_lu = feuille60.cell_value(k, 2)
    v_lu = feuille60.cell_value(k, 3)
    w_lu = feuille60.cell_value(k, 4)
    if u_lu != "NaN" and v_lu != "NaN" and w_lu != "NaN" :
        temps60.append(float(temps_lu)/1000) #temps en s
        u60.append(float(u_lu))
        v60.append(float(v_lu))
        w60.append(float(w_lu))
        norme60.append(np.sqrt(u60[-1]**2 + v60[-1]**2 + w60[-1]**2 ))
        if u60[-1]**2 > 0.001 : 
            direction60.append(np.arctan(v60[-1]/u60[-1])*180/np.pi)
        else :
            direction60.append(90)


if essai == 1:
    u30 = u30[Ndep:Ndep+20000]
    v30 = v30[Ndep:Ndep+20000]
    w30 = w30[Ndep:Ndep+20000]
    temps30 = temps30[Ndep:Ndep+20000]
    norme30 = norme30[Ndep:Ndep+20000]
    direction30 = direction30[Ndep:Ndep+20000]
    u45 = u45[Ndep:Ndep+20000]
    v45 = v45[Ndep:Ndep+20000]
    w45 = w45[Ndep:Ndep+20000]
    temps45 = temps45[Ndep:Ndep+20000]
    norme45 = norme45[Ndep:Ndep+20000]
    direction45 = direction45[Ndep:Ndep+20000]
    u60 = u60[Ndep:Ndep+20000]
    v60 = v60[Ndep:Ndep+20000]
    w60 = w60[Ndep:Ndep+20000]
    temps60 = temps60[Ndep:Ndep+20000]
    norme60 = norme60[Ndep:Ndep+20000]
    direction60 = direction60[Ndep:Ndep+20000]

# plt.plot(temps30,u30, color = "coral", linestyle = "--", label = "30 m, $u$")
# plt.plot(temps30,v30, color = "coral", linestyle = "-.", label = "30 m, $v$")
# plt.plot(temps30,w30, color = "coral", linestyle = ":", label = "30 m, $w$")
plt.plot(temps30,norme30, color = "coral", linestyle = "-", label = "30 m, $\\bar{U}$")
# plt.plot(temps45,u45, color = "limegreen", linestyle = "--", label = "45 m, $u$")
# plt.plot(temps45,v45, color = "limegreen", linestyle = "-.", label = "45 m, $v$")
# plt.plot(temps45,w45, color = "limegreen", linestyle = ":", label = "45 m, $w$")
plt.plot(temps45,norme45, color = "limegreen", linestyle = "-", label = "45 m, $\\bar{U}$")
# plt.plot(temps60,u60, color = "cyan", linestyle = "--", label = "60 m, $u$")
# plt.plot(temps60,v60, color = "cyan", linestyle = "-.", label = "60 m, $v$")
# plt.plot(temps60,w60, color = "cyan", linestyle = ":", label = "60 m, $w$")
plt.plot(temps60,norme60, color = "cyan", linestyle = "-", label = "60 m, $\\bar{U}$")
plt.xlabel("Temps (s)")
plt.ylabel("Vitesse (m/s)")
plt.title("Vitesses du vent enregistrées sur le mat")
plt.legend()

# plt.figure()
# plt.plot(temps30,direction30, color = "coral", linestyle = "-", label = "30 m, $\\phi_U$")
# plt.plot(temps45,direction45, color = "limegreen", linestyle = "-", label = "45 m, $\\phi_U$")
# plt.plot(temps60,direction60, color = "cyan", linestyle = "-", label = "60 m, $\\phi_U$")
# plt.xlabel("Temps (s)")
# plt.ylabel("Direction (°)")
# plt.title("Direction par rapport à l'est du vent")
# plt.legend()

#%% AFFICHAGE DES GRANDEURS UTILES
# Extraction des composantes stationnaires par lissage

fech = 10
Nstatio = 8*60*fech #la composante stationnaire moyenne est obtenue sur 8 minutes avec échantillage

def lissage(Lx, Ly, p):
    Lxout = []
    Lyout = []
    Lxout = Lx[p: -p]
    for index in range(p, len(Ly)-p):
        average = np.mean(Ly[index - p : index + p + 1])
        Lyout.append(average)
    return Lxout, Lyout

u30m = lissage(temps30, u30, Nstatio)[1]
v30m = lissage(temps30, v30, Nstatio)[1]
w30m = lissage(temps30, w30, Nstatio)[1]
norme30m = lissage(temps30, norme30, Nstatio)[1]
temps30m, direction30m = lissage(temps30, direction30, Nstatio)
u45m = lissage(temps45, u45, Nstatio)[1]
v45m = lissage(temps45, v45, Nstatio)[1]
w45 = lissage(temps45, w45, Nstatio)[1]
norme45m = lissage(temps45, norme45, Nstatio)[1]
temps45m, direction45m = lissage(temps45, direction45, Nstatio)
u60m = lissage(temps60, u60, Nstatio)[1]
v60m = lissage(temps60, v60, Nstatio)[1]
w60m = lissage(temps60, w60, Nstatio)[1]
norme60m = lissage(temps30, norme60, Nstatio)[1]
temps60m, direction60m = lissage(temps60, direction60, Nstatio)

plt.figure()
# plt.plot(temps30m,u30m, color = "coral", linestyle = "--", label = "30 m, $u$")
# plt.plot(temps30m,v30m, color = "coral", linestyle = "-.", label = "30 m, $v$")
# plt.plot(temps30mm,w30m, color = "coral", linestyle = ":", label = "30 m, $w$")
plt.plot(temps30m,norme30m, color = "coral", linestyle = "-", label = "30 m, $\\bar{U}$")
# plt.plot(temps45m,u45m, color = "limegreen", linestyle = "--", label = "45 m, $u$")
# plt.plot(temps45m,v45m, color = "limegreen", linestyle = "-.", label = "45 m, $v$")
# plt.plot(temps45m,w45m, color = "limegreen", linestyle = ":", label = "45 m, $w$")
plt.plot(temps45m,norme45m, color = "limegreen", linestyle = "-", label = "45 m, $\\bar{U}$")
# plt.plot(temps60m,u60m, color = "cyan", linestyle = "--", label = "60 m, $u$")
# plt.plot(temps60m,v60m, color = "cyan", linestyle = "-.", label = "60 m, $v$")
# plt.plot(temps60m,w60m, color = "cyan", linestyle = ":", label = "60 m, $w$")
plt.plot(temps60m,norme60m, color = "cyan", linestyle = "-", label = "60 m, $\\bar{U}$")
plt.xlabel("Temps (s)")
plt.ylabel("Vitesse (m/s)")
plt.title("Vitesses stationnaire du vent enregistrées sur le mat")
plt.legend()

plt.figure()
plt.plot(temps30m,direction30m, color = "coral", linestyle = "-", label = "30 m, $\\phi_U$")
plt.plot(temps45m,direction45m, color = "limegreen", linestyle = "-", label = "45 m, $\\phi_U$")
plt.plot(temps60m,direction60m, color = "cyan", linestyle = "-", label = "60 m, $\\phi_U$")
plt.xlabel("Temps (s)")
plt.ylabel("Direction (°)")
plt.title("Direction stationnaire par rapport à l'est du vent")
plt.legend()



#%% VERIFICATION DES LOIS LOGARITHMIQUES
# Détermination de la rugosité et du paramètre moyen

def prototype_log1(z, z0, U):
    return( U * np.log(z/z0) )

def prototype_log2(z, z0, U, zh):
    return( U * np.log((z-zh)/z0) )

def prototype_puis(z, z0, U, a):
    return( U * (z/z0)**a )

listz = [30, 45, 60]
vent1 = [min(norme30m), min(norme45m), min(norme60m)]
vent2 = [max(norme30m), max(norme45m), max(norme60m)]
vent3 = [np.mean(norme30m), np.mean(norme45m), np.mean(norme60m)]

initial_log1 = [1, 1]
initial_log2 = [1, 1, 1]
initial_puis = [1, 1, 1]

plt.figure()
plt.plot(listz, vent1, "bo")
plt.plot(listz, vent2, "go")
plt.plot(listz, vent3, "ro")
plt.xlabel("Hauteur $z$ (m)")
plt.ylabel("Vitesse moyenne du vent $\\bar{U}$ (m/s)")
plt.title("Essai des modèles explicites")

pfit, pcov = curve_fit(prototype_log1, listz, vent1, p0 = initial_log1)
print("Log 1 / Vent 1 (z0, U)",pfit)
plt.plot(listz, [prototype_log1(c, *pfit) for c in listz], "b--")
pfit, pcov = curve_fit(prototype_log1, listz, vent2, p0 = initial_log1)
print("Log 1 / Vent 2 (z0, U)",pfit)
plt.plot(listz, [prototype_log1(c, *pfit) for c in listz], "g--")
pfit, pcov = curve_fit(prototype_log1, listz, vent3, p0 = initial_log1)
print("Log 1 / Vent 3 (z0, U)",pfit)
plt.plot(listz, [prototype_log1(c, *pfit) for c in listz], "r--")

pfit, pcov = curve_fit(prototype_log2, listz, vent1, p0 = initial_log2)
print("Log 2 / Vent 1 (z0, U, zh)",pfit)
plt.plot(listz, [prototype_log2(c, *pfit) for c in listz], "b-.")
pfit, pcov = curve_fit(prototype_log2, listz, vent2, p0 = initial_log2)
print("Log 2 / Vent 2 (z0, U, zh)",pfit)
plt.plot(listz, [prototype_log2(c, *pfit) for c in listz], "g-.")
pfit, pcov = curve_fit(prototype_log2, listz, vent3, p0 = initial_log2)
print("Log 2 / Vent 3 (z0, U, zh)",pfit)
plt.plot(listz, [prototype_log2(c, *pfit) for c in listz], "r-.")

pfit, pcov = curve_fit(prototype_puis, listz, vent1, p0 = initial_puis)
print("Puis / Vent 1 (z0, U, a)",pfit)
plt.plot(listz, [prototype_puis(c, *pfit) for c in listz], "b:")
pfit, pcov = curve_fit(prototype_puis, listz, vent2, p0 = initial_puis)
print("Puis / Vent 2 (z0, U, a)",pfit)
plt.plot(listz, [prototype_puis(c, *pfit) for c in listz], "g:")
pfit, pcov = curve_fit(prototype_puis, listz, vent3, p0 = initial_puis)
print("Puis / Vent 3 (z0, U, a)",pfit)
plt.plot(listz, [prototype_puis(c, *pfit) for c in listz], "r:")


#%% DETERMINATION DES SPECTRES


direction_proj60 = np.mean(direction60m)*np.pi/180
norme_proj60 = np.mean(norme60m) # attention car prend en compte la partie verticale
ndata60 = len(temps60)
uproj60, uproj60f, vproj60f = [], [], []
for k in range(ndata60):
    uproj60.append(u60[k] * np.cos(direction_proj60) + v60[k] * np.sin(direction_proj60))
    uproj60f.append(uproj60[-1] - norme_proj60)
    vproj60f.append(-u60[k] * np.sin(direction_proj60) + v60[k] * np.cos(direction_proj60))

plt.figure()
plt.plot(temps60, uproj60f, label = "Vitesse $\\delta u$")
plt.plot(temps60, vproj60f, label = "Vitesse $\\delta v$")
plt.xlabel("Temps (s)")
plt.ylabel("Composante turbulente")
plt.legend()
plt.grid()


tfd = fft(w60)
N  =len(w60)
fe = 10
spectre = np.absolute(tfd)*2/N
freq = np.arange(N)*1.0/max(temps60)
plt.figure()
plt.plot(freq,spectre,'r')
plt.xlabel('f')
plt.ylabel('A')
plt.xscale("log")
plt.yscale("log")
plt.grid(True, which = "both", ls = "--", color= "0.7")
plt.axis([-0.1,fe/2,0,spectre.max()])
plt.grid()

# EN COURS
# https://www.f-legrand.fr/scidoc/docmml/numerique/tfd/spectre2/spectre2.html
# Il faudrait recalculer le spectre du vent qu'on propose dans le cube pour s'assurer qu'il suit bien le spectre !
