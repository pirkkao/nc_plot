#!/usr/bin/python3
#
# This is a simple python script to read in NetCDF data and produce
# 2D geographical maps of it with cartopy.
#
# Modify and play with this as you wish!
#

import mod_training as md
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Reload mod_training, needed when running in python-env
try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3
reload(md)



# ***********************************************************************
# DATA CONFIGURATION
# ***********************************************************************
#
# Create a dictionary for each data source you want to visualize
#
# 'exps'      : model forecast experiment identifier
#
# 'dates'     : forecast initialization dates
#
# 'types'     : what to plot
#                p000/ctrl : control forecast
#                p001      : ensemble member 1
#                p002      : ensemble member 2
#                ...       : ...
#                ensmean   : ensemble mean
#                ensstd    : ensemble standard deviation
#

main_dict=[]
main_dict.append({
        'exps'     : ["hands-on-6"],
        'dates'    : ["2017103100"],
        'types'    : ["ensmean","ensstd"]
        })


# ***********************************************************************
# PLOTTING CONFIGURATION
# ***********************************************************************
#
# Define variables to be plotted.
# 
# Fill in all the variables you want to plot into the pvars-array:
# [["MSL"], ["TP"]]
#
# For model fields on vertical levels, give the level number with the
# variable:
# [["T",1], ["T",10], ["MSL"]]
# 
# If only one variable is given, the dictionary will be filled with
# as many elements as there are data sources open.
#
# If plotting the cyclone tracks, all the elements will be filled
# with MSL according to the number of opened data sources.
#
pvars=[["W10M"]]



# Create a configuration dictionary
#
# 'fcsteps'    : forecast lengths to be plotted (choose "all"
#                for plotting all available timesteps)
#
# 'fig_name'   : unique identifier for the produced plot 
#                (.pdf added to this name)
#
# 'lonlat'     : area to be plotted 
#                "global" for global
#                "" for automatic area choosing according to data
#                [98.55,129.6,3.13796,23.75884] for an area
#
# 'minmax'     : find data min/max values over different forecast lengths
#                "abs" - get the absolute min/max values from all similar fields
#                "rel" - get the relative min/max for each individual field
#                "" - don't calculate min/max values
#
# 'plot_type'  : 
#                "2dmap" - plot a 2D map for each variable and forecast step defined
#                "mvar"  - plot multimple variables into the same plot,
#                          1st variable will be made with contour, 2nd with contourf
#                "track" - plot a simple tropical cyclone track
#

plot_dict={
    'fcsteps'    : range(14,26), # [10,11,12],
    'fig_name'   : "damrey_eps2",
    'lonlat'     : [104.,124,7.25,16.75],
    'minmax'     : "rel",
    'plot_type'  : "2dmap",

    # Define figure physical dimensions (size) and layout (ncol x nrow).
    # If left blank, default settings will used and ncol is defined to equal
    # number of variables to be plotted.
    'fig_size'   : (40,8),
    'fig_nrow'   : 1,
    'fig_ncol'   : 2,

    # Define number of contourf (cf_levs) and contour levels (c_levs), and
    # colour of contour lines.
    'fig_cf_levs': [],
    'fig_c_levs' : [],
    'fig_c_col'  : 'blue',

    # Associate a legend name for each variable. Automatically generated options:
    #    'type'
    #    'date'
    #    'exp'
    #    'datetype'
    #    'expdatetype'
    # or ["1","2","3","4"]
    'fig_legend' :'datetype',

    # Title, y- and x-labels
    'fig_title'  : [],
    'fig_ylabel' : [],
    'fig_xlabel' : [],

    # Change the cartopy projection
    'fig_proj'   : [],

    # Change observations used
    'fig_obs_track' : True, 
    'fig_obs_file'  : 'damrey_track.dat',
    'fig_obs_col'   : 'r',
    'fig_obs_buff'  : [],

    # Control plotting of additional map features
    'fig_features'  : True,
    }


# ***********************************************************************
# GET DATA
# ***********************************************************************

# Get variable list
plot_vars=md.create_vars(pvars,main_dict,plot_dict)

# Update plot_dict
plot_dict=md.configure_plot(plot_dict,plot_vars)

# Construct data paths
d_path = md.create_paths(main_dict,basepath="/wrk/ollinaho/DONOTREMOVE/public/")

# Fetch all data
dd = md.get_data_layer(d_path,parallel=False)

# Structure data for plotting
data_struct=md.structure_for_plotting(dd,plot_vars)

# Difference between two fields
#dd_st=[]
#dd_st.append(data_struct[0]-data_struct[1])
#data_struct=dd_st

# Get whole forecast len if requested
if plot_dict['fcsteps']=="all":
    plot_dict['fcsteps']=range(0,dd[0].sizes['time'])

# Get data min/max values from the forecast period
minmax=md.get_minmax_layer(plot_dict['minmax'],data_struct)


# ***********************************************************************
# PLOT
# ***********************************************************************

# Open a pdf to plot to
with PdfPages(plot_dict['fig_name']+'.pdf') as pdf:

    # PLOT 2D MAP
    if plot_dict['plot_type']=="2dmap":
        
        # Plot each forecast step into its own pdf page
        for fcstep in plot_dict['fcsteps']:

            # Call plotting
            #md.plot_data(fcstep,data_struct,lonlat=plot_dict['lonlat'],minmax=minmax,plot_vars=plot_vars)
            md.plot_data(fcstep,data_struct,plot_dict,plot_vars,minmax)

            # Save the plot to the pdf and open a new pdf page
            pdf.savefig()
            plt.close()


    # PLOT 2D MAP WITH MULTIPLE VARIABLES
    if plot_dict['plot_type']=="mvar":
        
        # Plot each forecast step into its own pdf page
        for fcstep in plot_dict['fcsteps']:

            # Call plotting
            md.plot_mvar(fcstep,data_struct,plot_dict,plot_vars,minmax)

            # Save the plot to the pdf and open a new pdf page
            pdf.savefig()
            plt.close()


    # PLOT TC TRACK
    elif plot_dict['plot_type']=="track":

        # Call plotting
        #md.plot_tctracks(data_struct,plot_dict['fcsteps'],plot_dict['lonlat'])
        md.plot_tctracks_and_pmin(data_struct,plot_dict,plot_vars)

        # Save the plot to the pdf and open a new pdf page
        pdf.savefig()
        plt.close()

