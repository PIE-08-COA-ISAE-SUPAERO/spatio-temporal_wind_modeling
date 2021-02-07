# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 17:39:18 2021

@author: 33658
"""


from vtk_functions import *
filename_wind = 'C:\\Users\\33658\\Desktop\\Cours 3A\\PIE\\data\\test_wind.vtk'
filename_surf =  'C:\\Users\\33658\\Desktop\\Cours 3A\\PIE\\data\\test_surf.vtk'

points = main(filename_wind,filename_surf)[0]
wind = main(filename_wind,filename_surf)[1]
surface = main(filename_wind,filename_surf)[2]
print ("3 first points:" , points[0:3])
print(" 3 first wind components:" ,wind[0:3] )
print("3 first surface points:", surface[0:3] )