import TurboMaxWindElitePro as wind

d = wind.wind()
d.import_wind_cube("G:/Mon Drive/PIE COA 08/Codes/PGM/spatio-temporal_wind_modeling/spatio-temporal_wind_modeling/data/Tests/exported_data_02-02-2021_13-15.json")

print(d.get_point(47,-113,2000))
#print(d.turbulence(47,-113,2000,10))
#d.profil_turbulence(47,-113,2000,1)
#d.plot_wind_cube([0,1000],[0,1000], [0,1000],True)
#d.plot_wind_surface("x",[47,-113],2000,True)
#d.plot_wind_surface("z",[47,-113], 2000, True)


