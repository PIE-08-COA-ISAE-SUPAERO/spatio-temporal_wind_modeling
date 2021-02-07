# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 17:38:06 2021

@author: 33658
"""

import numpy as np
import meshio



def main(filename_wind , filename_surf):
    """
    

    Parameters
    ----------
    filename_wind : String
        File containing the wind information, in a vtk format.
    filename_surf : String
        FIle containing the surface information, in a vtk format.

    Returns
    -------
    TYPE
        Returns 3 numpy arrays:
            The first containing 3 columns from the wind file, x y and z.
            The second one giving the 3 components of the wind at each position. (3 columns)
            The third are the ppints of the surface.
    """
    
    mesh = meshio.read(filename_wind)
    surface = meshio.read(filename_surf)
    return np.array(mesh.points),np.array(mesh.point_data['wind_vectors']),np.array(surface.points[:,2])
