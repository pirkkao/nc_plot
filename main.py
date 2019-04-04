#!/wrk/ollinaho/DONOTREMOVE/taito-conda-envs/plot2/bin/python3
#
# This is a simple python script to read in NetCDF data and produce
# 2D geographical maps of it with cartopy.
#
# Modify and play with this as you wish!
#
import mod_setup as msetup
import mod_data  as mdata
import mod_plot  as mplot
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

reload(msetup)
reload(mdata)
reload(mplot)
reload(mtrack)


# ***********************************************************************
# SETUP
# ***********************************************************************

exp="sv"

main_dict, pvars, operators, savescore, plot_dict2 = msetup.data_config(exp)

print("\nMAIN_DICT:","\n",main_dict)
print("\nPVARS:","\n",pvars)
print("\nOPERATORS:","\n",operators)
print("\nSCORE SAVING:","\n",savescore)
print("\nPLOTTING:","\n",plot_dict2,"\n")
#exit()

# ***********************************************************************
# AUTOMATIC CONFIGURATION
# ***********************************************************************

# Get plot dict
plot_dict = msetup.create_plot_dict()

# Get variable list
plot_vars = msetup.create_vars(pvars,main_dict,plot_dict)

# Update plot_dict
plot_dict = msetup.configure_plot(plot_dict,plot_vars)

# Include plot_dict2 TEMP SOLUTION?
plot_dict.update(plot_dict2)


# Construct data paths
d_path = msetup.create_paths(main_dict)

#exit()
# ***********************************************************************
# GET DATA
# ***********************************************************************

# Fetch data
data_struct = mdata.get_master(d_path,plot_vars,main_dict,operators,\
                               savescore,parallel=False)

# Get data min/max values from the forecast period
#minmax = mdata.get_minmax_layer(plot_dict['minmax'],data_struct)
minmax=[]

exit()
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

        # If CRPS divide CRPS and fair-CRPS plots
        if plot_dict['crps']:
            data= data_struct[0:-1:2]
            data2=data_struct[1:len(data_struct):2]
        else:
            data= data_struct
            data2=[]

        # Call plotting
        mplot.plot_scores3(plot_dict['fcsteps'],data,plot_dict,plot_vars,minmax)

        # Save the plot to the pdf and open a new pdf page
        pdf.savefig()
        plt.close()

        if data2:
            # Call plotting
            mplot.plot_scores3(plot_dict['fcsteps'],data2,plot_dict,plot_vars,minmax)

            # Save the plot to the pdf and open a new pdf page
            pdf.savefig()
            plt.close()


        # Call plotting
        #mplot.plot_scores2(plot_dict['fcsteps'],data_struct,plot_dict,plot_vars,minmax)

        # Save the plot to the pdf and open a new pdf page
        #pdf.savefig()
        #plt.close()
