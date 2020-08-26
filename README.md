# nc_plot
#

A python3 program to plot various aspects from NetCDF files 
using xarray. Main functionalities:
1) calculate skill scores from model output versus an analysis file
2) calculate time means of the skill scores
3) plot skill scores
4) plot 2D maps of various model output fields
5) plot maps of tropical cyclone tracks (simple MSLP minimum tracking)


Initialize from command line:
python3 main.py $name_of_config_file

e.g.
python3 main.py example_1A