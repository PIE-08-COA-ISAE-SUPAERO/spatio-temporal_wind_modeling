from typing import List


def main(input_path: str, simu_name: str):
    """
        Parameters
        ----------
        input_path : TYPE STRING
            Name of the input path, where you must put the .tif and  #Hour of collection of the data. Pleaqe type '00' for midnight, '06' for 6 AM and '22' for 10 PM.
        simu_name : TYPE STRING
            #DESCRIPTION. Altitude value for which we want to get the data. Please chose between 10m,20m,50m and 100m

        Returns
        -------
          : TYPE bool
            True, if it is works.
            False, if:
                The input folder doesn't exist, or
                The configuration file .json doesn't exist, or
                The terrain file .tif doesn't exist, or
                The wind file .??? doesn't exist, or
        """
    import os, json

    #Verify if the inputPath exist
    if os.path.exists(input_path):
        print("Input folder doesn't exist!")
        return False

    #Verify if the configuration file .json exist
    json_files: List[str] = [pos_json for pos_json in os.listdir(input_path) if pos_json.endswith('.json')]
    if json_files == []:
        print("Input .json file doesn't exist!")
        return False

    configFile: str = input_path + json_files[0]
    with open(configFile) as f:
         = json.load(f)


    #Verify if the terrain file .tif exist
    tif_files: List[str] = [pos_tif for pos_tif in os.listdir(input_path) if pos_tif.endswith('.tif')]
    if tif_files == []:
        print("Input .tif file doesn't exist!")
        return False
    else: tif_file: str = tif_files[0]

    #Verify if the wind file .??? exist
    wind_files: List[str] = [pos_wind for pos_wind in os.listdir(input_path) if pos_wind.endswith('.wind')]
    if wind_files == 0:
        print("Input .wind file doesn't exist!")
        return False
    else: wind_file: str = wind_files[0]


    # Create the simulation folder and move to him the input files# copie et cole les trois ficher de input dans le dossier simuName
    os.system("mkdir "+ "data/" + simu_name)
    os.system("move "+ tif_file + ", " + wind_file + "data/" + simu_name)

    ## lance windninja et assur que les output .vtk et .pdf vont être dans le dossier simuName




    configFile = ""

    call_WindNinja = "C:\\WindNinja\\WindNinja-3.5.3\\bin\\WindNinja_cli " + configFile + " --output_path " + outputFolder
    print(call_WindNinja)
    os.system(call_WindNinja)

    return True




    #simulationCase = input("Please enter the name of the case:\n")

    # flag = input("You have already defined the elevation file in the configuration file?(y/n)\n")
    # if flag == 'n':
    #     elevationFile = input("Please enter the name of the elevation file(with .tif):\n")
    #     elevationFile = "elevationFile = " + elevationFile
    #     outputFolder = "output_" + simulationCase + "_" + elevationFile
    # else:
    #     elevationFile = ""
    #    outputFolder = "output_" + simulationCase


    #
    # flag = input("You want to add any Simulation option?(y/n):\n")
    # newOption = ""
    # while flag == "y":
    #     newOption = input("Which option? Please write one of the ReadMe (name option argument)\n")
    #     newOption = " --" + newOption
    #     flag = input("You want to add any Simulation option?(y/n):\n")





