#!/usr/bin/python3
#
# This is a simple python script to read in NetCDF data and produce
# 2D geographical maps of it with cartopy.
#
# Modify and play with this as you wish!
#

import mod_data as mdata
import mod_plot as mplot
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Reload modules, needed when running in python-env
try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

reload(mdata)
reload(mplot)


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
#                ensmemb   : construct a list of ensemble member names,
#                            size indicated with 'nmem' : N.
# 

main_dict=[]
main_dict.append({
        'exps'     : ["hands-on-9"],
        'dates'    : ["2017110300"],
        'types'    : ["ensmean","p001","p002"],
        'nmem'     : 1,
        })

main_dict=mdata.update_main_dict(main_dict)


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
pvars=[["MSL"]]



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
    'fcsteps'    : range(5,15), # [10,11,12],(10,21)(20,31)(30,41)
    'fig_name'   : "testi",
    'lonlat'     : [106.,113,6.,16.],#[99.,129.,3.2,23.6][106.,113,6.,16.],
    'minmax'     : "rel",
    'plot_type'  : "2dmap",

    # Define figure physical dimensions (size) and layout (ncol x nrow).
    # If left blank, default settings will used and ncol is defined to equal
    # number of variables to be plotted.
    'fig_size'   : (12,14),
    'fig_nrow'   : [],
    'fig_ncol'   : [],

    # Define number of contourf (cf_levs) and contour levels (c_levs), and
    # colour of contour lines.
    'fig_cf_levs': 7,
    'fig_c_levs' : 30,
    'fig_c_col'  : 'magenta',

    # Track options. 
    #
    # 'fig_ens_predef' : Not in use
    # 'fig_ens_show'   : Show ensemble member tracks with solid lines
    # 'fig_ens_col'    : Use the same colour (defined here) for all ens members
    # 'fig_ens_buff'   : Halo size around ens member tracks
    # 'fig_ens_alpha'  : Alpha (transparency) of ens member halos
    # 'fig_ctrl_col'   : Colour for control member
    # 'fig_ensm_col'   : Colour for ensemble mean

    'fig_ens_predef' : False,
    'fig_ens_show'   : True,
    'fig_ens_col'    : [],
    'fig_ens_buff'   : [],
    'fig_ens_alpha'  : [],
    'fig_ctrl_col'   : [],
    'fig_ensm_col'   : 'magenta',

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
    'fig_obs_track'     : False, 
    'fig_obs_file'      : 'damrey_track.dat',
    'fig_obs_col'       : 'r',
    'fig_obs_buff'      : [],
    'fig_obs_match_time': True,

    # Control plotting of additional map features
    'fig_features'  : True,
    }


# ***********************************************************************
# GET DATA
# ***********************************************************************

# Get variable list
plot_vars=mdata.create_vars(pvars,main_dict,plot_dict)

# Update plot_dict
plot_dict=mdata.configure_plot(plot_dict,plot_vars)

# Construct data paths
d_path = mdata.create_paths(main_dict,basepath="/wrk/ollinaho/DONOTREMOVE/public/")

# Fetch all data
dd = mdata.get_data_layer(d_path,parallel=False)

# Structure data for plotting
data_struct=mdata.structure_for_plotting(dd,plot_vars)

# Deallocate data 
dd=[]

# Difference between two fields
#dd_st=[]
#dd_st.append(data_struct[0]-data_struct[1])
#data_struct=dd_st

# Get whole forecast len if requested
if plot_dict['fcsteps']=="all":
    plot_dict['fcsteps']=range(0,dd[0].sizes['time'])

# Get data min/max values from the forecast period
minmax=mdata.get_minmax_layer(plot_dict['minmax'],data_struct)


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
            #mplot.plot_data(fcstep,data_struct,lonlat=plot_dict['lonlat'],minmax=minmax,plot_vars=plot_vars)
            mplot.plot_data(fcstep,data_struct,plot_dict,plot_vars,minmax)

            # Save the plot to the pdf and open a new pdf page
            pdf.savefig()
            plt.close()


    # PLOT 2D MAP WITH MULTIPLE VARIABLES
    if plot_dict['plot_type']=="mvar":
        
        # Plot each forecast step into its own pdf page
        for fcstep in plot_dict['fcsteps']:

            # Call plotting
            mplot.plot_mvar(fcstep,data_struct,plot_dict,plot_vars,minmax)

            # Save the plot to the pdf and open a new pdf page
            pdf.savefig()
            plt.close()


    # PLOT TC TRACK
    elif plot_dict['plot_type']=="track":

        # Call plotting
        #mplot.plot_tctracks(data_struct,plot_dict['fcsteps'],plot_dict['lonlat'])
        mplot.plot_tctracks_and_pmin(data_struct,plot_dict,plot_vars)

        # Save the plot to the pdf and open a new pdf page
        pdf.savefig()
        plt.close()

