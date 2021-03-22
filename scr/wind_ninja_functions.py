# -*- coding: utf-8 -*-
"""Wind Ninja simulation.
This module contains all the functions required to launch the wind ninja simulation.
@author: A. Liberal Cavalcanti (March 2021)
"""
from typing import List

def main(input_path: str, simu_name: str, output_path: str):
    """Launch the wind ninja simulation and create the simulation files

        Parameters
        ----------
        input_path : String
            Name of the input path, where you must put the .tif and  #Hour of collection of the data. Pleaqe type '00' for midnight, '06' for 6 AM and '22' for 10 PM.
        output_path : String
            The path of the folder where every files is going to be stored
        simu_name : String
            The name of the simulation
        Returns
        -------
            bool
                True, if it is works.
                False, if:
                    The input folder doesn't exist, or
                    The configuration file .json doesn't exist, or
                    The terrain file .tif doesn't exist, or
                    The wind file .nc doesn't exist, or
        """
    import os, json, shutil

    # Verify if the inputPath exist
    if not os.path.exists(input_path):
        print("Input folder doesn't exist!")
        return False

    # Create the output folder
    output_path = output_path
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # # Create the simulation folder and move to him the input files
    # # copie et cole les trois ficher de input dans le dossier simu_name
    # data_path = './data'
    # simu_path = data_path + "/" + simu_name
    # if not os.path.exists(simu_path):
    #     os.mkdir("data/" + simu_name)

    # Verify if the configuration file .json exist
    json_files: List[str] = [pos_json for pos_json in os.listdir(input_path) if pos_json.endswith('.json')]
    # json_files_: List[str] = [pos_json for pos_json in os.listdir(simu_path) if pos_json.endswith('.json')]
    if not json_files:
        print("Input .json file doesn't exist!")
        return False

    # if not json_files_:
    json_file: str = json_files[0]
    # else:
    #     json_file: str = json_files_[0]
    
    json_path: str = output_path + simu_name + '.json'
    shutil.copyfile(input_path + json_file, json_path)
    
    with open(json_path) as f:
        config_text = json.load(f)
        version: str = config_text['def']['version']
        date: str = config_text['def']['date']
        mnt_file: bool = config_text['def']['mntFile']
        wind_file: bool = config_text['def']['windFile']
        wind_ninja_simulation = config_text['windNinjaSimulations']


    #Once the .json file is read it will get the mnt file's name and copy it in the simulation folder
    if mnt_file and input_path != output_path:
        tif_files: List[str] = [pos_tif for pos_tif in os.listdir(input_path) if pos_tif.endswith('.tif')]
        if not tif_files:
            print("Input .tif file doesn't exist!")
            return False
        else:
            tif_file: str = tif_files[0]
            shutil.copyfile(input_path + tif_file, output_path + tif_file)

    #Once the .json file is read it will get the wind file's and copy it in the simulation folder
    if wind_file and input_path != output_path:
        # Verify if the wind file .nc exist
        nc_files: List[str] = [pos_nc for pos_nc in os.listdir(input_path) if pos_nc.endswith('.nc')]
        if not nc_files:
            print("Input .nc file doesn't exist!")
            return False
        else:
            nc_file: str = nc_files[0]
            shutil.copyfile(input_path + nc_file, output_path + nc_file)

    # Write the .cfg, input config file of windNinja
    cfg_path = output_path + simu_name + '.cfg'
    with open(cfg_path, "w") as cfg_file:
        for param_key, param_values in wind_ninja_simulation.items():
            if param_key == 'fetch_elevation':
                cfg_file.write(param_key + '=' + str(output_path+param_values) + '\n')
            elif param_key != 'version' :
                cfg_file.write(param_key + '=' + str(param_values) + '\n')

    # Call the WindNinja_cli and call him in prompt command ligne
    # The outputs will be on the output_path
    windNinja_path = "C:/WindNinja/WindNinja-"+ version +"/bin/WindNinja_cli"

    try :
        windNinja_command: str = windNinja_path + " " + cfg_path + " --output_path " + output_path
        print(windNinja_command)
        os.system(windNinja_command)
    except:
        return False
    else :
        return True
