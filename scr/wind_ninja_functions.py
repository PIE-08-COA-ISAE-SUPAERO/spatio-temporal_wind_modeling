# -*- coding: utf-8 -*-
"""Wind Ninja simulation.
This module contains all the functions required to launch the wind ninja simulation.
@author: A. Cavalcanti Liberal (2020)
"""
from typing import List


def main(input_path: str, simu_name: str):
    """Launch the wind ninja simulation and create the simulation files

        Parameters
        ----------
        input_path : String
            Name of the input path, where you must put the .tif and  #Hour of collection of the data. Pleaqe type '00' for midnight, '06' for 6 AM and '22' for 10 PM.
        simu_name : String
            The path of the folder where every files is going to be stored
        Returns
        -------
            bool
                True, if it is works.
                False, if:
                    The input folder doesn't exist, or
                    The configuration file .json doesn't exist, or
                    The terrain file .tif doesn't exist, or
                    The wind file .??? doesn't exist, or
        """
    import os, json, shutil

    # Verify if the inputPath exist
    if not os.path.exists(input_path):
        print("Input folder doesn't exist!")
        return False

    # Create the simulation folder and move to him the input files
    # copie et cole les trois ficher de input dans le dossier simu_name
    data_path = './data'
    simu_path = data_path + "/" + simu_name
    if not os.path.exists(simu_path):
        os.mkdir("data/" + simu_name)

    # Verify if the configuration file .json exist
    json_files: List[str] = [pos_json for pos_json in os.listdir(input_path) if pos_json.endswith('.json')]
    json_files_: List[str] = [pos_json for pos_json in os.listdir(simu_path) if pos_json.endswith('.json')]
    if not json_files:
        print("Input .json file doesn't exist!")
        return False

    if not json_files_:
        json_file: str = json_files[0]
    else:
        json_file: str = json_files_[0]

    shutil.copyfile(input_path + json_file, simu_path + '/' + json_file)
    json_path: str = simu_path + '/' + json_file
    with open(json_path) as f:
        config_text = json.load(f)

        mnt_file = config_text['location']['mntFile']
        grib_file = config_text['location']['gribFile']
        for simu_id in range(len(config_text['windNinjaSimulations'])):
            chosen_simu = config_text['windNinjaSimulations'][simu_id]['name']
            if chosen_simu == simu_name:
                wind_ninja_simulation = config_text['windNinjaSimulations'][simu_id]
        config_data = [mnt_file, wind_ninja_simulation]

    print(config_data)
    print(wind_ninja_simulation.keys())

    # Verify if the terrain file .tif exist
    if mnt_file is True:
        tif_files: List[str] = [pos_tif for pos_tif in os.listdir(input_path) if pos_tif.endswith('.tif')]
        if not tif_files:
            print("Input .tif file doesn't exist!")
            return False
        else:
            tif_file: str = tif_files[0]
            shutil.copyfile(input_path + tif_file,
                            simu_path + '/' + tif_file)  # os.system("cp " + tif_file + ' ' + simu_path)

    if grib_file is True:
        # Verify if the wind file .??? exist
        wind_files: List[str] = [pos_wind for pos_wind in os.listdir(input_path) if pos_wind.endswith('.wind')]
        if not wind_files:
            print("Input .wind file doesn't exist!")
            return False
        else:
            wind_file: str = wind_files[0]
            shutil.copyfile(input_path + wind_file,
                            simu_path + '/' + wind_file)  # os.system("cp " + wind_file + ' ' + simu_path)

    # Write the .cfg
    cfg_name = simu_path + '/cli_' + simu_name + '.cfg'
    cfg_file = open(file=cfg_name, mode='w+')
    for param_key, param_values in wind_ninja_simulation.items():
        cfg_file.write(param_key + '=' + str(param_values) + '\n')

    # Create the output folder
    os.mkdir(simu_path + '/output')
    output_path = simu_path + '/output'

    ## lance windninja et assur que les output .vtk et .pdf vont être dans le dossier simuName

    call_WindNinja = "C:\\WindNinja\\WindNinja-3.5.3\\bin\\WindNinja_cli " + cfg_name + " --output_path " + output_path
    print(call_WindNinja)
    os.system(call_WindNinja)

    return True


