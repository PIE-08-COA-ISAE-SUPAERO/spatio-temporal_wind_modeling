# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 13:27:31 2021

@author: 33658
"""

# -*- coding: utf-8 -*-
"""PIE.ipynb

DONE BY GUEDIRA OUSSAMA

"""



from datetime import timedelta
import pandas as pd 
import math
import datetime
import pyproj
import pygrib



import xarray as xr

import numpy as np



























def speed(x, y):
    """
    

    Parameters
    ----------
    x : FLOAT
        
    y : FLOAT
        

    Returns root of the square of x plus the square of y. 
    Will be used to calculate the wind speed from u and v components
    -------
    TYPE
        FLOAT.

    """
    return math.sqrt(x**2 + y**2)



















def getDateMessage( grb):
    """
    

    Parameters
    ----------
    grb : grib message
        DESCRIPTION.

    Returns
    -------
    TYPE Datetime
        Given a message from a grib file, the function returns a datetime representing when the data has been collected.

    """
    
    try:
        return datetime.datetime(grb["year"],
                                           grb["month"], 
                                           grb["day"],
                                           grb["hour"],
                                           grb["minute"],
                                           grb["second"]) + timedelta(hours=grb["forecastTime"])
    except:
        print("Message does not exist, try another message with a lower number than grbs.messages")





    
















def chosen_param_man(file,altitude = 10,wind_component = "u-component of wind" ):
    """
    

    Parameters
    ----------
   file: TYPE STRING
   
    altitude : TYPE INT
        DESCRIPTION. Altitude value for which we want to get the data. Default value is 10.
    wind_component : TYPE String
        DESCRIPTION. Whether we want to get the u component of v component. The default is "u-component of wind".

    Returns
    -------
    TYPE grib message
        Returns the grib message (if exists, else returns 0) corresponding to the parameter and altitude chosen.

    """
    grbs = file
    try:
        grb = grbs.select(parameterName = wind_component,level = altitude)[0]
        return grb
    except:
        return 0






'''
print((grb.latitudes).shape)
print((grb.longitudes).shape)

print((grb.values.data.reshape(len(grb.latitudes),)).max())
print(grb.codedValues[0:5])'''












def grib_to_nc(file , path_name_output):
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
    

    
    
    













def main(hour_chosen,altitude_chosen,latitude_min,
         latitude_max, longitude_min,longitude_max, file, save_path):
    

    """
    

    Parameters
    ----------
    
    hour_chosen : String
        Hour of collection of the data. Pleaqe type '00' for midnight, '06' for 6 AM and '22' for 10 PM.
    altitude_chosen : TYPE INT
        DESCRIPTION. Altitude value for which we want to get the data. Please chose between 10m,20m,50m and 100m

    latitude_max : TYPE FLOAT
        maximum latitude value of points to keep
    latitude_min : TYPE FLOAT
        minimum latitude value of points to keep
    longitude_max : TYPE FLOAT
        maximum longitude value of points to keep
    longitude_min : TYPE FLOAT
        minimum longitude value of points to keep
    file : TYPE String
        Path and name of the grib file to use.
    save_path: TYPE STRING
    Path and Name to be given to the nc file that will be generated. For example type  
    'C:\\Users\\33658\\Desktop\\Cours_3A\\PIE\\data\\test_output.nc'


    Returns
    -------
    df : pandas dataframe
        Dataframe containing wind information for 
        every point having coordinates in the defined grid.
    The function generates an nc file depending on the entered parameters.


    """
    

    grib_to_nc(file, save_path)
    file = pygrib.open(file)

    
    grb = chosen_param_man(file,altitude = altitude_chosen,wind_component = "u-component of wind" )
    grb_2 = chosen_param_man(file,altitude = altitude_chosen,wind_component = "v-component of wind" )
    lst_v_component = list((grb_2.values.data.reshape(len(grb_2.latitudes),)))

    lst_lats = list(grb.latitudes)
    lst_longs = list(grb.longitudes)
    lst_u_component = list((grb.values.data.reshape(len(grb.latitudes),)))
    df = pd.DataFrame(list(zip(lst_lats, lst_longs , lst_u_component,lst_v_component)), 
               columns =['Latitudes', 'Longitudes' , "u-component of wind (m/s)" ,"v-component of wind (m/s)"]) 

    
    df = df[df["v-component of wind (m/s)"] != 9999]
    df = df[df["u-component of wind (m/s)"] != 9999] 
    df = df[df['Latitudes'].between(latitude_min - 0.001, latitude_max +0.001)]
    df = df[df['Longitudes'].between(longitude_min - 0.001, longitude_max + 0.001)]
    assert df.shape[0] > 0 , "No data for those values of latitude_max,latitude_min, longitude_max and longitude_min"


    df['wind speed (m/s)'] = df.apply(lambda x: speed(x['v-component of wind (m/s)'], x['u-component of wind (m/s)']), axis=1)
    df[["altitude (m)"]] = altitude_chosen
    df[["date"]] = getDateMessage(grb)

    return df

    

#main("03",20,45,45,10,11,'C:\\Users\\33658\\Downloads\\W_fr-meteofrance,MODEL,AROME+001+HP1+03H_C_LFPW_202103170000--.grib2', 'C:\\Users\\33658\\Desktop\\Cours_3A\\PIE\\data\\test_output.nc')


