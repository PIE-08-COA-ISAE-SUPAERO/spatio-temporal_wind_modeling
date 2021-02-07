# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 17:23:03 2021

@author: PIE COA 08:
     Sébastien Laigret
     Loïc Roux
     Alexandre Cavalcanti
     Oussama Guedira
     Quentin Abeille
"""

import numpy as np
# import tkinter as tk #Pour l'IHM
import json 

R_terre = 6371 #Rayon en km de la terre

class wind:
     
     def __init__(self):
          """
          The constructor of the class

          Returns
          -------
          None.

          """
          #Our cube with the position, wind speed and altitude of the surface
          self._wind_cube = {"Position" : [],
                             "Wind_speed" : [],
                             "Surface_altitude" : []}
          
          #The list of the point decomposed alongside every axis 
          self._list_point = {"x" : [],
                              "y" : [],
                              "z" : []}
          
          #The number of point inside the cube
          self._nb_points = 0
          
          #The height between 2 surfaces on the extrapolation part (outside of WindNinja first result)
          self._dz_extrapolation = 0
          
          #The highest point of the surface
          self._surface_elevation_max = 0
          
          #The date of the simulation
          self._date = ""
          
          #The location of the wind cube : the southwesternmost point of the cube
          self._location = (0,0)
          
          #The name of the file
          self._file_name = ""

     
     def _convert2dico(self):
          """
          Convert the obhect into a dictionnary 

          Returns
          -------
          wind_cube : dict
               The object as a dictionnary.

          """

          wind_cube = {
               "Date" : self._date,
               "Location" : self._location(),
               "Surface_elevation_max" : self._surface_elevation_max,
               "dz_extrapolation" : self._dz_extrapolation,
               "list_point" : self._list_point,
               "wind_cube" : self._wind_cube,
               "file_name" : self._file_name
               }
     
          return wind_cube
     
     def create_wind_cube(self, config_file):
          pass
     
     def import_wind_cube(self, file_name):
          pass
     
     def export_wind_cube(self):
          pass
     
     def get_data(self):
          return self._convert2dico()
     
     def get_point(self, latitude, longitude, altitude):
          """
          Return the wind parameters at the specific position requested, with interpolation if needed

          Parameters
          ----------
          latitude : double
               The latitude of the point wanted.
          longitude : double
               The longitude of the point wanted.
          altitude : double
               The altitude of the point wanted.

          Returns
          -------
          Tuple:
               u : double 
                    The west to east component of the wind speed
               v : double 
                    The south to north component of the wind speed
               w : double 
                    The vertical wind speed
               wind_speed : double 
                    The total wind speed
               wind_speed_flat: double 
                    The wind speed of the projection on the horizontal plan at the desired altitude
               direction : double 
                    The direction the wind came from, in ° 

          """
                    
          #Check if the point exists
          x_wanted, y_wanted = flat_distance_point((latitude, longitude), self._location)
          z_wanted = altitude
          
          min_z = min(self._list_point["z"])
          max_z = max(self._list_point["z"])
          
          if self._list_point["x"][-1] < x_wanted or self._list_point["y"][-1] < y_wanted or z_wanted < min_z or max_z < z_wanted :
               print("This point does not belong in the cube")
               return 0
          
          #Get the smallest cube possible 
          x_min, x_max = smallest_interval(x_wanted, self._list_point["x"])
          y_min, y_max = smallest_interval(y_wanted, self._list_point["y"])
          z_min, z_max = smallest_interval(z_wanted, self._list_point["z"])
          
          #Get the distance for the interpolation
          dx = x_max - x_min
          dy = y_max - y_min
          dz = z_max - z_max
          
          #Get the 2 extremeties : lowest southwesternmost and highest northesternmost
          cube_min = (x_min, y_min, z_min)
          cube_max = (x_max, y_max, z_max)
          
          #Get the index in the list of those two
          i_min = self._wind_cube["Position"].index(cube_min)
          i_max = self._wind_cube["Position"].index(cube_max)
          
          #Get the associated wind speeds
          u_min, v_min, w_min = self._wind_cube["Wind_speed"][i_min]
          u_max, v_max, w_max = self._wind_cube["Wind_speed"][i_max]
     
          #Create the interpolation
          u = u_min + (x_wanted - x_min) * (u_max - u_min) / dx
          v = v_min + (y_wanted - y_min) * (v_max - v_min) / dy
          w = w_min + (z_wanted - z_min) * (w_max - w_min) / dz
          
          #Compute the wind speeds and direction
          wind_speed = np.sqrt(u**2 + v**2 + w**2)
          wind_speed_flat = np.sqrt(u**2 + v**2)
          
          direction = np.rad2deg(np.arccos(v / wind_speed_flat))
          #Correction to bring if the wind comes from the west instead of the east
          if u < 0 : direction *= -1
          
          return u, v, w, wind_speed, wind_speed_flat, direction
     
     
     
     
def smallest_interval(x, array):
     """
     Get the smallest interval possible which contains the value x in the array 
     
     Parameters
     ----------
     x : float
          The value who want to find
     array : n-array
          The array which contains the value to look into.
     
     Returns
     -------
     Tuple or 0 if the interval was not found
          The smallest interval that contains the value x
     
     """
     n = len(array)
     
     for i in range(n-1):
          if array[i] <= x and x <= array[i+1]:
               return (array[i], array[i+1])
     
     return 0 
def flat_distance_point(point1, point2):
     """
     Compute the distance as a straight line between 2 geographical coordinates

     Parameters
     ----------
     point1 : tuple
          The latitude and longitude of the first point.
     point2 : tuple
          The latitude and longitude of the second point.

     Returns
     -------
     dx : double
          The distance over the x axis (along the east-west direction).
     dy : double
          The distance over the y axis (along the north-south direction).

     """
     
     (lat1, long1) = point1
     (lat2, long2) = point2
     
     dlat = abs(lat1-lat2)/2
     dlong = abs(long1-long2)/2
     
     dx = 2 * R_terre * np.sin(dlong)
     dy = 2 * R_terre * np.sin(dlat)
     
     return dx, dy