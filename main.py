#!/usr/bin/python3
#
# This is a simple python script to read in NetCDF data and produce
# 2D geographical maps of it with cartopy.
#
# Modify and play with this as you wish!
#

import mod_data as mdata
import mod_plot as mplot
import mod_tctrack as mtrack
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
reload(mtrack)


# ***********************************************************************
# DATA CONFIGURATION
# ***********************************************************************

exp="eda+sv"

main_dict, pvars, operators = mdata.data_config(exp)

main_dict = mdata.update_main_dict(main_dict)

print("MAIN_DICT:")
print(main_dict)
print()

print("PVARS:")
print(pvars)
print()

print("OPERATORS:")
print(operators)
print()

#exit()

# ***********************************************************************
# PLOTTING CONFIGURATION
# ***********************************************************************

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
    'fcsteps'    : range(5,9), # [10,11,12],(10,21)(20,31)(30,41)
    'fig_name'   : "testi",
    'lonlat'     : "global", #[106.,113,6.,16.],#[99.,129.,3.2,23.6][106.,113,6.,16.],
    'minmax'     : "rel",
    'plot_type'  : "score",

    # Define figure physical dimensions (size) and layout (ncol x nrow).
    # If left blank, default settings will used and ncol is defined to equal
    # number of variables to be plotted.
    'fig_size'   : (12,14),
    'fig_nrow'   : 2,
    'fig_ncol'   : 1,

    # Define number of contourf (cf_levs) and contour levels (c_levs), and
    # colour of contour lines.
    'fig_cf_levs': 20,
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
    'fig_legend' :'type',

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
plot_vars = mdata.create_vars(pvars,main_dict,plot_dict)

# Update plot_dict
plot_dict = mdata.configure_plot(plot_dict,plot_vars)

# Construct data paths
d_path = mdata.create_paths(main_dict)

getdata=False
score_fname='crps_fair_T850_Nx.nc'

if getdata:
    # Fetch all data
    data_struct = mdata.get_data_layer(d_path,plot_vars,parallel=False)

    # Data operations
    data_struct, nam_list = mdata.structure_for_plotting2(data_struct,main_dict,operators)

    mdata.save_score_data(score_fname,data_struct,nam_list)

else:
    data_struct=mdata.get_score_data(score_fname)

#print(data_struct)
exit()

# Get whole forecast len if requested
if plot_dict['fcsteps']=="all":
    plot_dict['fcsteps']=range(0,dd[0].sizes['time'])

# Get data min/max values from the forecast period
#minmax = mdata.get_minmax_layer(plot_dict['minmax'],data_struct)
minmax=[]

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
        mtrack.plot_tctracks_and_pmin(data_struct,plot_dict,plot_vars)

        # Save the plot to the pdf and open a new pdf page
        pdf.savefig()
        plt.close()


    # PLOT SCORES
    elif plot_dict['plot_type']=="score":

        # Call plotting
        #mplot.plot_scores_crps_vs_fair(plot_dict['fcsteps'],data_struct,plot_dict,plot_vars,minmax)

        # Save the plot to the pdf and open a new pdf page
        #pdf.savefig()
        #plt.close()


        # Call plotting
        mplot.plot_scores3(plot_dict['fcsteps'],data_struct[0:10:2],plot_dict,plot_vars,minmax)

        # Save the plot to the pdf and open a new pdf page
        pdf.savefig()
        plt.close()


        # Call plotting
        mplot.plot_scores3(plot_dict['fcsteps'],data_struct[1:10:2],plot_dict,plot_vars,minmax)

        # Save the plot to the pdf and open a new pdf page
        pdf.savefig()
        plt.close()


        # Call plotting
        mplot.plot_scores2(plot_dict['fcsteps'],data_struct,plot_dict,plot_vars,minmax)

        # Save the plot to the pdf and open a new pdf page
        pdf.savefig()
        plt.close()
