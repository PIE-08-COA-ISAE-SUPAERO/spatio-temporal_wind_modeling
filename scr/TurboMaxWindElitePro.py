# -*- coding: utf-8 -*-
"""Wind class
This module contains the wind class and every function linked to it
@author: PIE COA 08 (2020)
     Sébastien Laigret
     Loïc Roux
     Alexandre Cavalcanti
     Oussama Guedira
     Quentin Abeille
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd

import json
import os 

#import wind_ninja_functions as wn_function
import vtk_functions as vtk
import extrapolation 
import plots

#%% Global variable
R_terre = 6371 #Earth radius (km)

#Some variable for turbulence
A_Z = 0.177
B_Z = 0.00823
C_Z = 1.2

PENTE_PRIN = 5/6
D_PRIN = 70.7
PENTE_LATW = 11/6
D_LATW = 188.4

TIME_MAX = 900

CORREC_GLOBAL = 1/3
CORREC_ALTI = 0.1
RATIO_PRIN = 0.11
RATIO_LAT = 0.11
RATIO_VERT = 0.06

#%% The class
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
          
          #The number of surface layer on the extrapolation part (outside of WindNinja first result)
          self._nb_layer_extrapolation = 6
          
          #The highest point (m) you want the simulation to go
          self._elevation_max = 0
          
          #The date of the simulation
          self._date = ""
          
          #The location of the wind cube : the southwesternmost point of the cube
          self._location = [0,0]
          
          #The name of the folder
          self._folder_name = ""

     def __str__(self):
          """Overloading of the print function
          Returns
          -------
          String
              The string to print
          """          
          txt = ""
          for key, val in self._convert2dico().items():
               if key == "Location":
                    txt += key + " : (" + str(val[0]) + ", " + str(val[1]) +"),\n"
               elif key == "list_point":
                    txt += key  + " : " + "{x,y,z}, " + str(self._nb_points) + " elts,\n"
               elif key == "wind_cube":
                    txt += key  + " : " + "{position, wind speeds, surface altitude}, " + str(self._nb_points) + " elts,\n"
               else:
                    txt += key + " : " + str(val) + ",\n"

          return txt


     def _convert2dico(self):
          """
          Convert the obhect into a dictionnary 
          Returns
          -------
          wind_cube : dict
               The object as a dictionnary.
          """
          #We need to convert into array the np-array
          wind_cube_temp = {}
          for key, array in self._wind_cube.items():
               wind_cube_temp[key] = array.tolist()

          list_point_temp = {}
          for key, array in self._list_point.items():
               list_point_temp[key] = array.tolist()

          wind_cube = {
               "Date" : self._date,
               "Location" : self._location,
               "Elevation_max" : self._elevation_max,
               "nb_layer_extrapolation" : self._nb_layer_extrapolation,
               "nb_points": self._nb_points,
               "folder_name" : self._folder_name, 
               "list_point" : list_point_temp,
               "wind_cube" : wind_cube_temp
               }

          return wind_cube
     
     def create_wind_cube(self, input_folder_path, simulation_folder_name):
          """Create the wind cube by launching the simulation
          Parameters
          ----------
          input_folder_path : String
              The folder where the input file (config_file : .json, .tif, and the wind file are stored)
          simulation_folder_name : String
              The name of the folder where the data is going to be stored in the DATA folder
          Returns
          -------
               flag :Bool
                    A boolean which testify if the simulation has worked 
          """   

          flag = False

          #We get the config file to extract some information         
          #Check if there is only 1 file
          files = file_list_by_extension(input_folder_path, ".json")          
          if len(files) > 1 :
               print("There is more than 1 .json file, please only keep one")
               return flag

          config_file = files[0]
          with open(input_folder_path + "/" + config_file, "r") as read_file:
               data = json.load(read_file)
          
          #We store the information we have
          self._nb_layer_extrapolation = data["extrapolation"]["nb_layer_extrapolation"]
          self._elevation_max = data["extrapolation"]["elevation_max"]
          self._location = (data["location"]["latitude"],data["location"]["longitude"])
          self._folder_name = simulation_folder_name
          self._date = data["period"]["start"]

          #We launch the Wind Ninja simulation
          assert wn_function.main(input_folder_path, simulation_folder_name)

          #We get the .vtk file
          #Check if there is only 1 file
          files = file_list_by_extension(simulation_folder_name, ".vtk")
          if len(files) != 2 :
               print("There is an error in the simulation")
               return flag
          
          wind_file  = surf_file = simulation_folder_name + '/'

          if "_surf" in files[0] :     
               wind_file += files[1]
               surf_file += files[0]
          else:
               wind_file += files[0]
               surf_file += files[1]

          points, wind, surface = vtk.main(wind_file, surf_file)
          extrap_field = extrapolation.main(points, wind, surface, self._elevation_max, self._nb_layer_extrapolation)

          #We construct the dictionnary cube
          self._wind_cube["Position"] = extrap_field[: , 0:3]
          self._wind_cube["Wind_speed"] = extrap_field[: , 3:6]
          self._wind_cube["Surface_altitude"] = extrap_field[: , 6]

          self._list_point["x"] = np.unique(self._wind_cube["Position"][:,0])
          self._list_point["y"] = np.unique(self._wind_cube["Position"][:,1])
          self._list_point["z"] = np.unique(self._wind_cube["Position"][:,2])

          self._nb_points = len(extrap_field)

          #We export the cube
          assert self.export_wind_cube()

          #The simulation has worked 
          print("The wind cube is created")
          flag = True
          return flag 

     def import_wind_cube(self, file_name):
          """Import the data from a save file 
          Parameters
          ----------
          file_name : String
              The complete path of the save file
          """          
          with open(file_name, 'r') as read_file :
               data = json.load(read_file)
          
          self._date = data["Date"]
          self._location = data["Location"]
          self._elevation_max = data["Elevation_max"]
          self._nb_layer_extrapolation = data["nb_layer_extrapolation"]
          self._folder_name = data["folder_name"]
          self._nb_points = data["nb_points"]

          #We need to convert the list into np-array
          self._wind_cube["Position"] = np.array(data["wind_cube"]["Position"])
          self._wind_cube["Wind_speed"]= np.array(data["wind_cube"]["Wind_speed"])
          self._wind_cube["Surface_altitude"]= np.array(data["wind_cube"]["Surface_altitude"])

          self._list_point["x"] = np.array(data["list_point"]["x"])
          self._list_point["y"] = np.array(data["list_point"]["y"])
          self._list_point["z"] = np.array(data["list_point"]["z"])
     
     def export_wind_cube(self):
          """
          Create the .json file and stock it 
          Returns
          -------
          int
               1 if the file has correctly been exported.
          """
          name = self._folder_name + "/exported_data_" + self._date + ".json"
          with open(name, "w+") as f:
               json.dump(self._convert2dico(), f)
               print("Export done")
          return True
               
     def get_data(self):
          """Return the dictionnary of the object
          Returns
          -------
          dict
              The dictionary which represents the object
          """          
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
          #Get the distance the cube uses
          x_wanted, y_wanted = flat_distance_point((latitude, longitude), self._location)
          z_wanted = altitude
          
          #Check if the point is in the x and y boundaries
          if self._list_point["x"][-1] < x_wanted or self._list_point["y"][-1] < y_wanted :
               print("This point does not belong in the cube")
               return 0

          #Get the surface altitude as a function of x and y
          map_surface = [ [ self._wind_cube["Position"][i, 0], self._wind_cube["Position"][i, 1], self._wind_cube["Surface_altitude"][i] ] for i in range(self._nb_points)]
          map_surface = np.array(map_surface)
          
          #Get the smallest cube (x,y) possible 
          x_min, x_max = smallest_interval(x_wanted, self._list_point["x"])
          y_min, y_max = smallest_interval(y_wanted, self._list_point["y"])
          
          #Find the surface in this small cube
          surface_alt = extrapolation.get_Zlist_pos(x_max, y_max, map_surface)[1]
          surface_alt = np.unique(surface_alt)[0]

          #Check the z boundary
          min_z = surface_alt
          max_z = min_z + self._elevation_max          
          if z_wanted < min_z:
               print("This point is below ground : " + str(z_wanted) + "m VS " + str(min_z) + "m")
               return 0

          if max_z < z_wanted :
               print("This point is to high above the ground : " + str(z_wanted) + "m VS " + str(max_z) + "m")
               return 0
          
          z_min, z_max = smallest_interval(z_wanted, self._list_point["z"])


          #Get the distance for the interpolation
          dx = x_max - x_min
          dy = y_max - y_min
          dz = z_max - z_min
          
          #Get the 2 extremeties : lowest southwesternmost and highest northesternmost
          cube_min = [x_min, y_min, z_min]
          cube_max = [x_max, y_max, z_max]
          
          #Get the index in the list of those two
          i_min = np_array_index(cube_min, self._wind_cube["Position"])
          i_max = np_array_index(cube_max, self._wind_cube["Position"])
          
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
          
          return u, v, w, wind_speed, wind_speed_flat, direction_plan(u, v)

     def turbulence(self, latitude, longitude, altitude, time):
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
          
          u, v, w, _, wind_speed_flat, direction = self.get_point(latitude, longitude, altitude)
          
          #We get the altitude of the selected geographical point
          map_surface = [ [ self._wind_cube["Position"][i, 0], self._wind_cube["Position"][i, 1], self._wind_cube["Surface_altitude"][i] ] for i in range(self._nb_points)]
          map_surface = np.array(map_surface)
          x, y = flat_distance_point((latitude, longitude), self._location)

          #Get the smallest cube (x,y) possible 
          _, x_max = smallest_interval(x, self._list_point["x"])
          _, y_max = smallest_interval(y, self._list_point["y"])

          surface_alt = extrapolation.get_Zlist_pos(x_max, y_max, map_surface)[1]
          surface_alt = np.unique(surface_alt)[0]

          _, _, _, _, wind_10, _ = self.get_point(latitude, longitude, surface_alt + 10)

          fs = 10 # frequence sampling
          N = int(time * fs) # Max number of modes
          
          # Table of random Fourier coefficients
          altitude_correction = max((CORREC_ALTI - 1)/(1000 - 0)*altitude + 1, CORREC_ALTI)
          s_base = wind_10 * altitude_correction * CORREC_GLOBAL
          s_prin = RATIO_PRIN * s_base
          s_late = RATIO_LAT * s_base
          s_w = RATIO_VERT * s_base 
          
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
          direction_prin = 270 - direction_plan(wind_prin, wind_late)
          direction_base = 270 - direction
                                   
          direction_prov = direction_base + direction_prin
          u = magnitude_plan * np.cos(direction_prov * np.pi/180)
          v = magnitude_plan * np.sin(direction_prov * np.pi/180)
          w = wind_w
          direction = direction_plan(u, v)  
          
          # Plot of the windrose
          fig = plt.figure()
          ax = fig.add_subplot(111, projection="polar")
          title = "Windrose (" + str(round(latitude,2)) + "°," + str(round(longitude,2)) + "°) : " + str(round(magnitude_plan,1)) + " m/s, " + str(round(direction,1)) + "°"
          ax.set_title(title)
          ax.set_rmin(0)
          ax.set_rmax(max(25,magnitude_plan + 5))
          ax.set_thetalim(0, 2*np.pi)
          ax.set_theta_direction(-1)
          ax.set_theta_zero_location("N")
          plt.arrow(direction*np.pi/180, magnitude_plan + 0.5, 0, -magnitude_plan, \
                    width = 0.1, lw = 1, length_includes_head = True, head_width=1, \
                         head_length=3, shape = "full", overhang = 0.5)
          plt.show()
          return(u, v, w, magnitude, magnitude_plan, direction)


     def profil_turbulence(self, latitude, longitude, altitude, timestep):
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
          
          u, v, w, _, wind_speed_flat, direction = self.get_point(latitude, longitude, altitude)
          
          #We get the altitude of the selected geographical point
          map_surface = [ [ self._wind_cube["Position"][i, 0], self._wind_cube["Position"][i, 1], self._wind_cube["Surface_altitude"][i] ] for i in range(self._nb_points)]
          map_surface = np.array(map_surface)
          x, y = flat_distance_point((latitude, longitude), self._location)

          #Get the smallest cube (x,y) possible 
          _, x_max = smallest_interval(x, self._list_point["x"])
          _, y_max = smallest_interval(y, self._list_point["y"])

          surface_alt = extrapolation.get_Zlist_pos(x_max, y_max, map_surface)[1]
          surface_alt = np.unique(surface_alt)[0]

          _, _, _, _, wind_10, _ = self.get_point(latitude, longitude, surface_alt + 10)
          
          fs = 1/timestep # frequence sampling
          TEMPS_MAX = 900 # observation of the profile during 15 minutes
          N = int(TEMPS_MAX * fs) # Max number of modes
          
          # Table of random Fourier coefficients
          # Table of random Fourier coefficients
          X_prin = []
          X_late = []
          X_w = []
            
          # Ponderation of the impact of turbulence thanks to global coefficients
          altitude_correction = max((CORREC_ALTI - 1)/(1000 - 0)*altitude + 1, CORREC_ALTI)
          s_base = wind_10 * altitude_correction * CORREC_GLOBAL
          s_prin = RATIO_PRIN * s_base
          s_late = RATIO_LAT * s_base
          s_w = RATIO_VERT * s_base 
            
            
          for k in range(N):
              frequency = k/(2*N) * fs
              s_prink = np.sqrt(TIME_MAX/(2*np.pi) * s_prin**2 * spectre_prin(frequency,\
                                                        wind_speed_flat, altitude))
              s_latek = np.sqrt(TIME_MAX/(2*np.pi) * s_late**2 * spectre_late(frequency,\
                                                            wind_speed_flat, altitude))
              s_wk = np.sqrt(TIME_MAX/(2*np.pi) * s_w**2 * spectre_w(frequency, \
                                                            wind_speed_flat, altitude))
              x_prin = rd.normal(0, s_prink)
              x_late = rd.normal(0, s_latek)
              x_w = rd.normal(0, s_wk)
              X_prin.append(x_prin)
              X_late.append(x_late)
              X_w.append(x_w)
            
          # Computing the values of the wind in the time by reverse FFT
          list_time = [k * timestep for k in range(N)]
          list_prin = []
          list_late = []
          list_w = []    
          list_u = []
          list_v = []
            
          list_prin2 = np.fft.ifft(X_prin)
          list_late2 = np.fft.ifft(X_late)
          list_w2 = np.fft.ifft(X_w)
        
          for k in range(N) :
                
              wind_prin = wind_speed_flat + 2*np.pi * fs * list_prin2[k].real
              wind_late = 0 + 2*np.pi * fs * list_late2[k].real
              wind_w = w + 2*np.pi * fs * list_w2[k].real
                
              list_prin.append(wind_prin)
              list_late.append(wind_late)
              list_w.append(wind_w)
                
              magnitude_plan = np.sqrt(wind_prin**2 + wind_late**2)      
                
              direction_prin = 270 - direction_plan(wind_prin, wind_late)
              direction_base = 270 - direction
              direction_new = direction_base + direction_prin
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
          plt.show()
          plt.figure()
          plt.plot(list_time, list_u, label = "West/East wind $U_x$")
          plt.plot(list_time, list_v, label = "South/North wind $U_y$")
          plt.plot(list_time, list_w, label = "Vertical wind $U_z$")
          plt.xlabel("Time (s)")
          plt.ylabel("Speed of the wind (m/s)")
          plt.title("Evolution in time of the wind speed")
          plt.grid()
          plt.legend()
          plt.show()
          
          # Output as the data computed for prediction
          return(list_time, list_prin, list_late, list_u, list_v, list_w)

     def plot_wind_surface(self, axis, coord, alt, plot):
          """
          Calculates a wind profile or a wind surface from the wind_cube.
          Parameters
          ----------
          axis : string
               Type of interpoltion to do. "z" for a wind surface, "x" or "y" for a
               wind profile.
          coord : Narry of float
               Coordinates of the point for the wind profile.
          alt : float
               Altitude for the wind surface.
          plot : Bool
               Activates the plot option.
          Returns
          -------
          X_mesh : Narray of floats
               3D-mesh for the interpolated x coordinates.
          Y_mesh : Narray of floats
               3D-mesh for the interpolated x coordinates.
          Z_mesh : Narray of floats
               3D-mesh for the interpolated x coordinates.
          Uinterp : Narray of floats
               3D-mesh for the interpolated wind speed component along x-axis.
          Vinterp : Narray of floats
               3D-mesh for the interpolated wind speed component along y-axis.
          Winterp : Narray of floats
               3D-mesh for the interpolated wind speed component along z-axis.
          Sinterp : Narray of floats
               3D-mesh for the interpolated surface altitude.
          """
          coordxy_tuple = flat_distance_point(self._location, coord)
          return plots.plot_wind_surface(self._wind_cube, axis, [coordxy_tuple[0],coordxy_tuple[1]], alt, plot)
     
     def plot_wind_cube(self, xlim, ylim, zlim, plot):
          """
          Calculates a wind_cube inside the wind_cube and plots it if needed.
          Parameters
          ----------
          xlim : Narray of floats
               Numpy array of size 2 containing the limits of the x-axis range.
          ylim : Narray of floats
               Numpy array of size 2 containing the limits of the y-axis range.
          zlim : Narray of floats
               Numpy array of size 2 containing the limits of the z-axis range.
          plot : Bool
               Activates the plot or not.
               
          Returns
          -------
          X_mesh : Narray of floats
               3D-mesh for the interpolated x coordinates.
          Y_mesh : Narray of floats
               3D-mesh for the interpolated x coordinates.
          Z_mesh : Narray of floats
               3D-mesh for the interpolated x coordinates.
          Uinterp : Narray of floats
               3D-mesh for the interpolated wind speed component along x-axis.
          Vinterp : Narray of floats
               3D-mesh for the interpolated wind speed component along y-axis.
          Winterp : Narray of floats
               3D-mesh for the interpolated wind speed component along z-axis.
          Sinterp : Narray of floats
               3D-mesh for the interpolated surface altitude.
          """
          return plots.plot_wind_cube(self._wind_cube, xlim, ylim, zlim, plot)

#---------------------------------------------------------
# #%% AUXILIARY
#---------------------------------------------------------
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
     
     if x <= array[0] : return (0, array[0])

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
    
    lu = altitude/((A_Z+B_Z*altitude)**(C_Z))
    S = 4*lu/speed * 1/((1 + D_PRIN*(frequency*lu/speed)**2)**(PENTE_PRIN))
    
    return(S)


def spectre_late(frequency, speed, altitude):
    
    """
    Return the DSP of the lateral wind speed
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
    
    lv = altitude/((A_Z+B_Z*altitude)**(C_Z))
    S = 4*lv/speed * (1 + D_LATW*(2*frequency*lv/speed)**2)/((1 + D_PRIN*(2*frequency*lv/speed)**2)**(PENTE_LATW))
    
    return(S)


def spectre_w(frequency, speed, altitude):
    
    """
    Return the DSP of the vertical wind speed
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
    
    lw = altitude
    S = 4*lw/speed * (1 + D_LATW*(2*frequency*lw/speed)**2)/((1 + D_PRIN*(2*frequency*lw/speed)**2)**(PENTE_LATW))
    
    return(S)


def file_list_by_extension(file_dir, file_extension):
     return  [_ for _ in os.listdir(file_dir) if _.endswith(file_extension)] 

def np_array_index(point, list_point):
     """Find the index of a np-array inside a 2D np-array
     Parameters
     ----------
     point : np-array
         the np-array to find 
     list_point : 2D np-array
         The 2D np-array where to find the point inside
     Returns
     -------
     Integer
         The index of the point inside list_point
     """     
     [x, y, z] = point
     
     cond_x = list_point[:,0] == x
     cond_y = list_point[:,1] == y
     cond_z = list_point[:,2] == z

     cond = np.multiply(cond_x, cond_y)  
     cond = np.multiply(cond, cond_z)
     
     if len(np.argwhere(cond==True)) == 0:
          return 0
     else :    
          return np.argwhere(cond==True)[0,0]

def direction_plan(u, v):
    """
    Return the direction of the origin of the wind (where it comes from)
    Parameters
    ----------
    u : double 
        The west to east component of the wind speed
    v : double 
        The south to north component of the wind speed
    Returns
    -------
    direction : double 
         The angle, in degree, of the coming wind
    """
    
    magnitude = np.sqrt(u**2 + v**2)
    angle = 180/np.pi * np.arccos(abs(u)/magnitude)
    if u > 0 and v > 0 :
        direction = 270 - angle
    if u > 0 and v < 0 :
        direction = 270 + angle
    if u < 0 and v < 0 :
        direction = 90 - angle
    if u < 0 and v > 0 :
        direction = 90 + angle
    
    return(direction%360)