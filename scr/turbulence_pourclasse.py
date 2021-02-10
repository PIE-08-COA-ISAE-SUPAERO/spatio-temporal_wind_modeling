# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 11:54:39 2021

@author: Sébastien
VERSION NON DEFINITIVE POUR CLASSE :
    - j'enlèverai les exemples et rajouterai les self
    - non prise en compte de l'hypothèse de Taylor
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd
import time

#%% FONCTIONS AUXILIAIRES (spectre, directions principales)


def spectre_prin(frequency, speed, altitude):
    
    """
    Return the DSP of the principal wind speed
    Parameters
    ----------
    frequency : double
         The frequency of the component expected.
    speed : double
         The stationnary speed of the wind at the point wanted.
    altitude : double
         The altitude of the point wanted.  
    Returns
    -------
    S : double 
         The value of the DSP for these parameters
    """
    
    lu = altitude/((0.177+0.00823*altitude)**(1.2))
    S = 4*lu/speed * 1/((1 + 70.7*(frequency*lu/speed)**2)**(5/6))
    
    return(S)


def spectre_late(frequency, speed, altitude):
    
    lv = altitude/((0.177+0.00823*altitude)**(1.2))
    S = 4*lv/speed * (1 + 188.4*(2*frequency*lv/speed)**2)/((1 + 70.7*(2*frequency*lv/speed)**2)**(11/6))
    
    return(S)


def spectre_w(frequency, speed, altitude):
    
    lw = altitude
    S = 4*lw/speed * (1 + 188.4*(2*frequency*lw/speed)**2)/((1 + 70.7*(2*frequency*lw/speed)**2)**(11/6))
    
    return(S)



#%% FONCTION DE PROFIL TURBULENT

def turbulence(latitude, longitude, altitude, time):
    
    """
    Return the wind parameters at the specific position requested with a random turbulence component
    Parameters
    ----------
    latitude : double
         The latitude of the point wanted.
    longitude : double
         The longitude of the point wanted.
    altitude : double
         The altitude of the point wanted.
    time : double
         The relative time (seconds) since the initial time    
    Returns
    -------
    Tuple:
         u : double 
              The west to east component of the wind speed
         v : double 
              The south to north component of the wind speed
         w : double 
              The vertical wind speed
         magnitude : double 
              The total wind speed
         magnitude_plan : double 
              The wind speed of the projection on the horizontal plan at the desired altitude
         direction : double 
              The direction the wind came from, in ° 
    """
    
    # Attention, time est un délai relatif à la mesure stationnaire de Windninja, pas une heure absolue.
    #u, v, w, wind_speed, wind_speed_flat, direction = self.get_point(self, \
    #                                            latitude, longitude, altitude)
    # u, v, w, wind_speed, wind_speed_flat, direction = 10, 5, 1, np.sqrt(10**2 \
    #             + 5**2 + 1), np.sqrt(10**2 + 5**2), np.arctan(5/10)*180/np.pi
    u, v, w, wind_speed, wind_speed_flat, direction = 20, 0, 0, 20, 20, 0
    altitude = 70
    
    fs = 10 # frequence sampling
    N = int(time * fs) # Max number of modes
    
    # Table of random Fourier coefficients
    s_prin, s_late, s_w = 1, 1, 1 #
    wind_prin, wind_late, wind_w = wind_speed_flat, 0, w
    
    for k in range(N):
        
        frequency = k/(2*N) * fs
        s_prink = np.sqrt(time/(2*np.pi) * s_prin**2 * spectre_prin(frequency,\
                                                    wind_speed_flat, altitude))
        s_latek = np.sqrt(time/(2*np.pi) * s_late**2 * spectre_late(frequency,\
                                                    wind_speed_flat, altitude))
        s_wk = np.sqrt(time/(2*np.pi) * s_w**2 * spectre_w(frequency, \
                                                    wind_speed_flat, altitude))
        x_prin = rd.normal(0, s_prink)
        x_late = rd.normal(0, s_latek)
        x_w = rd.normal(0, s_wk)
        wind_prin = wind_prin + 2*np.pi*fs/N * x_prin * np.cos(2*np.pi*k/N*time)
        wind_late = wind_late + 2*np.pi*fs/N * x_late * np.cos(2*np.pi*k/N*time)
        wind_w = wind_w + 2*np.pi*fs/N * x_w * np.cos(2*np.pi*k/N*time)

    # Deviation with regards the principal axis
    magnitude = np.sqrt(wind_prin**2 + wind_late**2 + wind_w**2)
    magnitude_plan = np.sqrt(wind_prin**2 + wind_late**2)
    direction_prin = np.arctan(wind_late/wind_prin)*180/np.pi
                            
    direction = direction + direction_prin
    u = magnitude_plan * np.cos(direction * np.pi/180)
    v = magnitude_plan * np.sin(direction * np.pi/180)
    w = wind_w    
    
    return(u, v, w, magnitude, magnitude_plan, direction)


def profil_turbulence(latitude, longitude, altitude, timestep):
    
    """
    Return the evolution of the speed of the wind at a specific point during 15 minutes
    Parameters
    ----------
    latitude : double
         The latitude of the point wanted.
    longitude : double
         The longitude of the point wanted.
    altitude : double
         The altitude of the point wanted.
    timestep : double
         The time step (in seconds) between each value of the wind    
    Returns
    -------
    Tuple:
         list_time : list of doubles
              The time scale used for measures
         list_prin : list of doubles
              The evolution of the principal perturbated wind speed
         list_late : list of doubles
              The evolution of the lateral perturbated wind speed
         list_u : list of doubles
              The evolution of the west to east component of the wind speed
         list_v : list of doubles
              The evolution of the south to north component of the wind speed
         list_w : list of doubles
              The evolution of the vertical wind speed 
    """
    
    # Attention, time est un délai relatif à la mesure stationnaire de Windninja, pas une heure absolue.
    #u, v, w, wind_speed, wind_speed_flat, direction = self.get_point(self, \
    #                                            latitude, longitude, altitude)
    u, v, w = 15, 5, 2
    wind_speed, wind_speed_flat, direction = np.sqrt(u**2 + v**2 + w**2), np.sqrt(u**2 + v**2), np.arctan(v/u)*180/np.pi
    #u, v, w, wind_speed, wind_speed_flat, direction = 20, 0, 0, 20, 20, 0
    
    fs = 1/timestep # frequence sampling
    TEMPS_MAX = 900 # observation of the profile during 15 minutes
    N = int(TEMPS_MAX * fs) # Max number of modes
    
    # Table of random Fourier coefficients
    s_prin, s_late, s_w = 1, 1, 1 # to ponderate the impact of turbulence
    X_prin = []
    X_late = []
    X_w = []
    
    for k in range(N):
        frequency = k/(2*N) * fs
        s_prink = np.sqrt(TEMPS_MAX/(2*np.pi) * s_prin**2 * spectre_prin(frequency,\
                                                    wind_speed_flat, altitude))
        s_latek = np.sqrt(TEMPS_MAX/(2*np.pi) * s_late**2 * spectre_late(frequency,\
                                                    wind_speed_flat, altitude))
        s_wk = np.sqrt(TEMPS_MAX/(2*np.pi) * s_w**2 * spectre_w(frequency, \
                                                    wind_speed_flat, altitude))
        x_prin = rd.normal(0, s_prink)
        x_late = rd.normal(0, s_latek)
        x_w = rd.normal(0, s_wk)
        X_prin.append(x_prin)
        X_late.append(x_late)
        X_w.append(x_w)
    
    # Computing the values of the wind in the time
    list_time = []
    list_prin = []
    list_late = []
    list_u = []
    list_v = []
    list_w = []
    
    for k in range(int(TEMPS_MAX/timestep)):
        
        t = k * timestep
        wind_prin, wind_late, wind_w = wind_speed_flat, 0, w
        for kk in range(N) :
            wind_prin = wind_prin + 2*np.pi*fs/N * X_prin[kk] * np.cos(2*np.pi*kk/N*t)
            wind_late = wind_late + 2*np.pi*fs/N * X_late[kk] * np.cos(2*np.pi*kk/N*t)
            wind_w = wind_w + 2*np.pi*fs/N * X_w[kk] * np.cos(2*np.pi*kk/N*t)
        list_time.append(t)
        list_prin.append(wind_prin)
        list_late.append(wind_late)
        list_w.append(wind_w)
        magnitude_plan = np.sqrt(wind_prin**2 + wind_late**2)
        direction_prin = np.arctan(wind_late/wind_prin)*180/np.pi
        direction_new = direction + direction_prin
        u = magnitude_plan * np.cos(direction_new * np.pi/180)
        v = magnitude_plan * np.sin(direction_new * np.pi/180)
        list_u.append(u)
        list_v.append(v)
        
    # Output as a graph of wind evolution in time
    plt.figure()
    plt.plot(list_time, list_prin, label = "Principal wind $\\bar{U}+\\tilde{u}$")
    plt.plot(list_time, list_late, label = "Lateral wind $\\tilde{v}$")
    plt.plot(list_time, list_w, label = "Vertical wind $\\bar{W}+\\tilde{w}$")
    plt.xlabel("Time (s)")
    plt.ylabel("Speed of the wind (m/s)")
    plt.title("Evolution in time of the wind speed")
    plt.grid()
    plt.legend()
    plt.figure()
    plt.plot(list_time, list_u, label = "West/East wind $U_x$")
    plt.plot(list_time, list_v, label = "South/North wind $U_y$")
    plt.plot(list_time, list_w, label = "Vertical wind $U_z$")
    plt.xlabel("Time (s)")
    plt.ylabel("Speed of the wind (m/s)")
    plt.title("Evolution in time of the wind speed")
    plt.grid()
    plt.legend()
    
    # Output as the data computed for prediction
    return(list_time, list_prin, list_late, list_u, list_v, list_w)

t1 = time.time()
turbulence(1, 1, 70, 100)
t2 = time.time()
list_time, list_prin, list_late, list_u, list_v, list_w = profil_turbulence(1, 1, 70, 0.1)
t3 = time.time()
print(t2-t1, t3-t2)
print("u,v,w", np.mean(list_u), np.mean(list_v), np.mean(list_w), "U, vp", np.mean(list_prin), np.mean(list_late))