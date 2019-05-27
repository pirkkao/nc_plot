#!/wrk/ollinaho/DONOTREMOVE/taito-conda-envs/plot2/bin/python3
#
# This is a simple python script to read in NetCDF data and produce
# either skill scores against observations or make 2D geographical maps 
# of it with cartopy.
#
# Modify and play with this as you wish!
#
import mod_setup   as msetup
import mod_data    as mdata
import mod_plot    as mplot
import mod_tctrack as mtrack
import sys

# Reload modules, needed when running in python-env
try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

reload(msetup)
reload(mdata)
reload(mplot)
reload(mtrack)


# ***********************************************************************
# SETUP
# ***********************************************************************

# Experiment name
exp="ex1"

# Type of action (data/plot)
exptyp="plot"

# Directory from which to fetch config (leave blank for default)
expdir=""

# Overwrite if command line arguments given
if len(sys.argv)>1:
    exp=str(sys.argv[1])
if len(sys.argv)>2:
    exptyp=str(sys.argv[2])
if len(sys.argv)>3:
    expdir=str(sys.argv[3])

# Call setup
main_dict, pvars, operators, dataoper, savescore, plot_dict2 = msetup.data_config(exp,exptyp,expdir)


#exit()
# ***********************************************************************
# AUTOMATIC CONFIGURATION
# ***********************************************************************

# Get plot dict
plot_dict = msetup.create_plot_dict()
plot_dict.update(plot_dict2)

# Get variable list
plot_vars = msetup.create_vars(pvars,main_dict,plot_dict)

# Update plot_dict
plot_dict = msetup.configure_plot(plot_dict,plot_vars)

# Update plot_dict (again)
plot_dict.update(plot_dict2)

# Construct data paths
d_path = msetup.create_paths(main_dict,plot_vars)


#exit()
# ***********************************************************************
# GET DATA
# ***********************************************************************

# Fetch data
data_struct = mdata.get_master(d_path,plot_vars,main_dict,operators,\
                               dataoper,savescore,plot_dict,parallel=False)

# Get data min/max values from the forecast period
minmax = mdata.get_minmax_layer(plot_dict,data_struct)


#exit()
# ***********************************************************************
# PLOT
# ***********************************************************************

# Call plot
if plot_dict['plot_type'] != 'none':
    mplot.plot_master(data_struct,plot_dict,plot_vars,operators,minmax)
