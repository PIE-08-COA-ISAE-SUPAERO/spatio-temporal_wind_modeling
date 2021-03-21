# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 20:08:35 2021

@author: 33658
"""


import xarray as xr

def main(file , path_name_output):
    """
    

    Parameters
    ----------
    file : TYPE String
        Grib file to read.
    path_name_output : TYPE String
        Path and name of the file to provide

    Returns
    -------
    None.

    """
    ds = xr.open_dataset(file, engine='cfgrib')
    ds.to_netcdf(path_name_output)
    

    