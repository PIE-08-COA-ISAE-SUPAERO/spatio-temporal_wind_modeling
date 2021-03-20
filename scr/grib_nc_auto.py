# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 20:09:59 2021

@author: 33658
"""



import datetime
import webbrowser
import time
import glob
#import cfgrib


import xarray as xr






def get_grib_from_web(paquet_chosen,hour_chosen,downloads_path):
    """
    

    Parameters
    ----------
    paquet_chosen : String 
        Either 'SP1' for surface data or 'HP1' for higher altitudes. SP1 will give data at an altitude of 10 meters,
        HP1 will give data at altitudes of 20,50 and 100 meters.
        
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
    
    
    print('If the download did not start, please copy and paste the following link in another browser: ', site)

    while file_and_path not in glob.glob(downloads_path+'\\*.grib2'):
        print("Still downloading")
        time.sleep(5)
        print()
    print('file downloaded successfully')
    
    return 





def get_grib_from_path(paquet_chosen,hour_chosen, downloads_path):
    """
    

    Parameters
    ----------
    paquet_chosen : String 
        Either 'SP1' for surface data or 'HP1' for higher altitudes.
        
    hour_chosen : String
        Hour of collection of the data. Pleaqe type '00' for midnight, '06' for 6 AM and '22' for 10 PM.
    downloads_path : String
        Path where the grib files are stored. Should be the same as the path where the downloads are stored

    Returns
    -------
    file : String
        The function returns the filename of the grib file completing all the requirements. If the file does not exist, it returns False

    """
    files_list = glob.glob(downloads_path+'\\*.grib2')
    for file in files_list:
        if file.split('+')[2] == paquet_chosen:
            if file.split('_')[1][-3:-1] == hour_chosen:
                if file.split('_')[4][0:8] == datetime.datetime.now().strftime("%Y%m%d"):
                    return file
    return False


def main(paquet_chosen,hour_chosen, downloads_path, path_name_output):
    """
    

    Parameters
    ----------
    paquet_chosen : String 
        Either 'SP1' for surface data or 'HP1' for higher altitudes.
        
    hour_chosen : String
        Hour of collection of the data. Pleaqe type '00' for midnight, '06' for 6 AM and '22' for 10 PM.
    downloads_path : String
        Path where the grib files are stored. Should be the same as the path where the downloads are stored
    path_name_output : TYPE String
        Path and name of the file to provide

    Returns
    -------
    None.

    """
   
    
    get_grib_from_web(paquet_chosen,hour_chosen,downloads_path)
    file = get_grib_from_path(paquet_chosen,hour_chosen, downloads_path)
    
    ds = xr.open_dataset(file, engine='cfgrib')
    ds.to_netcdf(path_name_output)
    return 
    





main('HP1',"03",'C:\\Users\\33658\\Downloads', 'C:\\Users\\33658\\Desktop\\Cours_3A\\PIE\\data\\test_output.nc')

                