# -*- coding: utf-8 -*-
"""PIE.ipynb

DONE BY GUEDIRA OUSSAMA

"""



from datetime import timedelta
import pandas as pd 
import math
import datetime
import webbrowser
import time
import glob
import pyproj
import pygrib



import xarray as xr

import numpy as np








def open_file_grib(file_name):
    """
    

    Parameters
    ----------
    file_name : String
        Name of a grib file (file format: grib2).

    Returns
    -------
    bool
        Returns true if the file is not readable (does not exist), else returns False.

    """
    try:
        pygrib.open(file_name)
        print("File downloaded successfully")
        return False
    except:
        print("Still downloading")
        return True









def get_grib_from_web(paquet_chosen,hour_chosen,downloads_path):
    """
    

    Parameters
    ----------
    paquet_chosen : String 
        Either 'SP1' for surface data or 'HP1' for higher altitudes.
        
    hour_chosen : String
        Hour of collection of the data. Pleaqe type '00' for midnight, '06' for 6 AM and '22' for 10 PM.
    downloads_path : String
        Default path of where the downloads are stored.
        
        
    The function will directly download the required file and check if the operation has been successfully completed.
    Will stop only if it is the case.

    Returns
    -------
    None.

    """
    now = datetime.datetime.now()
    paquet = paquet_chosen
    hour = hour_chosen
    date = now.strftime("%Y-%m-%d")
    date_format_2 = now.strftime("%Y%m%d")
    
    
    
    site = "https://donneespubliques.meteofrance.fr/?fond=donnee_libre&token=__5yLVTdr-sGeHoPitnFc7TZ6MhBcJxuSsoZp6y0leVHU__&model=AROME&format=grib2&grid=0.01&grid2=0.01&package={}&time={}H&referencetime={}T00%3A00%3A00Z".format(paquet,hour,date)
    webbrowser.open(site)
    file_name = 'W_fr-meteofrance,MODEL,AROME+001+{}+{}H_C_LFPW_{}0000--.grib2'.format(paquet_chosen,hour_chosen,date_format_2)
    file_and_path = downloads_path + '\\' + file_name
    while open_file_grib(file_and_path) :
        time.sleep(5)
        
    return 
















def get_grib_from_path(paquet_chosen,hour_chosen,path):
    """
    

    Parameters
    ----------
    paquet_chosen : String 
        Either 'SP1' for surface data or 'HP1' for higher altitudes.
        
    hour_chosen : String
        Hour of collection of the data. Pleaqe type '00' for midnight, '06' for 6 AM and '22' for 10 PM.
    path : String
        Path where the grib files are stored. Should be the same as the path where the downloads are stored

    Returns
    -------
    file : String
        The function returns the filename of the grib file completing all the requirements.

    """
    files_list = glob.glob(path+'\\*.grib2')
    for file in files_list:
        if file.split('+')[2] == paquet_chosen:
            if file.split('_')[1][-3:-1] == hour_chosen:
                if file.split('_')[4][0:8] == datetime.datetime.now().strftime("%Y%m%d"):
                    return file
                













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





    
















def chosen_param(HP1,SP1,altitude = 10,wind_component = "u-component of wind" ):
    """
    

    Parameters
    ----------
    HP1 : TYPE grib file read with pygrib.open
        First grib file giving data at altitudes of 20,50 and 100 meters.
    SP1 : TYPE grib file read with pygrib.open
        Second grib file giving data at altitude of 10 meters. 
        !!!Please make sure that both files contain data for the same datetime!!!
    altitude : TYPE INT
        DESCRIPTION. Altitude value for which we want to get the data. Default value is 10.
    wind_component : TYPE String
        DESCRIPTION. Whether we want to get the u component of v component. The default is "u-component of wind".

    Returns
    -------
    TYPE grib message
        Returns the grib message (if exists, else returns 0) corresponding to the parameter and altitude chosen.

    """
    if altitude == 10:
        grbs = SP1
    elif altitude in[20,50,100]:
        grbs = HP1
    else:
        print("Please chose altitude of 10,20,50 or 100m")
        return 0
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


























def main(hour_chosen,altitude_chosen,latitude_min,
         latitude_max, longitude_min,longitude_max, downloads_path, save_path):
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
    downloads_path : TYPE String
        Path where the grib files are stored. Should be the same as the path where the downloads are stored
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
    assert altitude_chosen in [10,20,50,100] , "Please choose altitude between 10m,20m,50m and 100m."
    assert hour_chosen in ["0"+str(i) for i in range(0,10)]+ [str(i) for i in range(10,25)] , "Please type '00' for    midnight, '06' for 6 AM and '22' for 10 PM." 
    
    get_grib_from_web("HP1",hour_chosen,downloads_path)
    get_grib_from_web("SP1",hour_chosen,downloads_path)
    
    HP1_file = pygrib.open(get_grib_from_path('HP1',hour_chosen,downloads_path))
    SP1_file = pygrib.open(get_grib_from_path('SP1',hour_chosen,downloads_path))


    
    grb = chosen_param(HP1_file,SP1_file,altitude = altitude_chosen,wind_component = "u-component of wind" )
    grb_2 = chosen_param(HP1_file,SP1_file,altitude = altitude_chosen,wind_component = "v-component of wind" )
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
    net = df[[ 'Latitudes', 'Longitudes' , "u-component of wind (m/s)" ,"v-component of wind (m/s)"]].to_xarray()
    new_filename_2 = save_path
    print ('saving to ', new_filename_2)
    net.to_netcdf(path=new_filename_2)
    dataset_1 = xr.open_dataset(new_filename_2)
    
    for key in dataset_1.keys():
        print ('checking %s ' % key)
        print ('-- identical in dataset_2 and net_dataset : %s'\
               % np.allclose(dataset_1[key], net[key], equal_nan=True))
    return df

    

#main("03",10,45,45,10,11,'C:\\Users\\33658\\Downloads', 'C:\\Users\\33658\\Desktop\\Cours_3A\\PIE\\data\\test_output.nc')


