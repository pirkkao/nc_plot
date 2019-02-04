
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
import seaborn as sns
import shapely.geometry as sgeom
from shapely.geometry import Point
from datetime import datetime,timedelta
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D as Line

# import from local
from mod_util import plot_legend, plot_features

# Don't display FutureWarnings
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")



def plot_tctracks_and_pmin(data_struct,plot_dict,plot_vars):
    "Plot TC tracks for given data"

    fcsteps=plot_dict['fcsteps']
    lonlat=plot_dict['lonlat']

    # Create a figure
    fig=plt.figure(figsize=(12,15))

    ax1=fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax2=fig.add_subplot(212)


    # Get date of first data file for later
    dt0 = datetime.strptime(plot_vars[0]['dates'], "%Y%m%d%H")     

    # Create colours
    cols=col_list_data(data_struct,plot_dict,plot_vars)

    
    ifc_diff=0
    icol=0
    for data in data_struct:
        # Plot forecast low track
        pmins=create_tc_track(ax1,data,plot_vars[icol],fcsteps,cols[icol],buff=plot_dict['fig_ens_buff'],\
                                  ens_show=plot_dict['fig_ens_show'],buff_alpha=plot_dict['fig_ens_alpha'])

        # Construct x-axis for pmin plot based on forecast initialization date
        xax=[x*3 for x in fcsteps]

        # Check date of currently opened data
        dt = datetime.strptime(plot_vars[icol]['dates'], "%Y%m%d%H")        

        if dt != dt0:
            # Get the time difference in time steps [3h]
            ifc_diff=int((dt-dt0).seconds/3600./3 +(dt-dt0).days*8)

            xax=[(x+ifc_diff)*3 for x in fcsteps]
            print("Calculated difference in timesteps for data source number ",str(icol)," is",ifc_diff)

        ax2.plot(xax,pmins,color=cols[icol])

        icol+=1


    # Plot Damrey track
    if plot_dict['fig_obs_match_time']:
        xax_obs,obs=tc_plot(ax1,'damrey_track.dat','red',buff=plot_dict['fig_obs_buff'],\
                                match_date_to=[dt0,fcsteps])
    else:
        tc_plot(ax1,'damrey_track.dat','red',buff=plot_dict['fig_obs_buff'],\
                    match_date_to=[])

    # Plot Damrey observed pressure
    if plot_dict['fig_obs_match_time']:
        ax2.plot(xax_obs,obs,color='r')
    

    # Change x-tick properties
    ax2.set_xticks([x*3 for x in range(fcsteps[0],fcsteps[-1]+ifc_diff,2)])
    ax2.set_xlabel("Hours from "+str(dt0))

    # Legend
    plot_legend(cols,[],plot_vars=plot_vars,plot_dict=plot_dict)

    # Plot map features
    plot_features(ax1,country=True,lam=True,lonlat=lonlat)





def create_tc_track(ax,data,plot_var,fcsteps,icol,buff=[],ens_show=False,buff_alpha=[]):
    "Find and construct the TC track from the data"

    lons=[]
    lats=[]
    pmin=[]

    for fcstep in fcsteps:
        # Do a slice north of 10N since there is another low in the area
        ilon,ilat,ipmin=find_pressure_min(data.isel(time=fcstep).sel(lat=slice(20.,10.)))

        if not ilon==[]:
            if len(ilon)>1:
                print("More than one minimum ",ilon)
            lons.append(float(ilon[0]))
            lats.append(float(ilat[0]))
            pmin.append(float(ipmin))
        else:
            pmin.append(float('nan'))


    # Calculate distance between points
    #trackpoint_prev=[]
    #tmp_lons=[]
    #tmp_lats=[]
    #ilon=
    #for trackpoint in zip(lons,lats):
    #   if trackpoint_prev:
    #        distance = trackpoint_prev.distance(Point(trackpoint))
    #        
    #        if distance < 4.:
    #            tmp_lons.append([])
    #            tmp_lats.append([])
    #
    #    trackpoint_prev=Point(trackpoint)
    #
    #    ilon+=1

    # Turn the lons and lats into a shapely LineString
    track = sgeom.LineString(zip(lons, lats))

    # Plot the track
    if buff and plot_var['ens']=="member":
        alpha=0.2
    else:
        alpha=0.7

    if plot_var['ens']=="ctrl":
        ax.add_geometries([track], ccrs.PlateCarree(),
                          facecolor='none',edgecolor='k',linewidth=3.0,alpha=alpha)

    elif plot_var['ens']!="member" or (plot_var['ens']=="member" and ens_show):
        ax.add_geometries([track], ccrs.PlateCarree(),
                          facecolor='none',edgecolor=icol,linewidth=3.0,alpha=alpha)

    # Plot a buffer around the track
    if buff and plot_var['ens']=="member":
        track_buffer=track.buffer(buff)

        if not buff_alpha:
            buff_alpha=0.15

        ax.add_geometries([track_buffer], ccrs.PlateCarree(),
                          facecolor=icol, alpha=buff_alpha)

    return pmin



def find_pressure_min(data):
    "Simply find the coordinates of the lowest pressure reading"

    # Find the value
    pmin=data.min().values

    # Filter values with over threshold pressure
    if pmin < 100550.:
        # Find the data point
        ploc=data.where(data==pmin,drop=True)

        return ploc['lon'].values,ploc['lat'].values,pmin

    else:
        return [],[],[]




def tc_plot(ax,data_path,edgecolor,buff,match_date_to=[]):
    "Plot TC track with given lat-lon lists"

    # Get data 
    lats,lons = tc_track(data_path)
    dates,pmins = tc_depth(data_path)

    if match_date_to:
        dt0=match_date_to[0]

        d_start=dt0+timedelta(hours=match_date_to[1][0]*3)
        d_stop =dt0+timedelta(hours=match_date_to[1][-1]*3)

        print("\nDamrey requested track first data point:", d_start)
        print("                       last  data point:", d_stop,"\n")

        # Obs first date
        dt1=dates[0]

        # Get the time difference in time steps [6h]
        ifc_diff=int((d_start-dt1).seconds/3600./6 +(d_start-dt1).days*4)
        ifc_last=int((d_stop -dt1).seconds/3600./6 +(d_stop -dt1).days*4)

        # Also diff to initialization date
        ifc_init=int((dt1-dt0).seconds/3600./6 +(dt1-dt0).days*4)

        print("I think I need to adjust obs time by ",ifc_diff, \
                  " elements and end at element ", ifc_last,"\n")

        print("New observation time span is: ",dates[ifc_diff])
        print("                              ",dates[ifc_last],"\n")

        print(ifc_init)

        lats=lats[ifc_diff:ifc_last]
        lons=lons[ifc_diff:ifc_last]


        xax=[x*6 for x in range(ifc_diff+ifc_init,ifc_last+ifc_init+1)]
        print(xax)
        print(pmins[ifc_diff:ifc_last+1])
        

    # Turn the lons and lats into a shapely LineString
    track = sgeom.LineString(zip(lons, lats))


    # Plot the track
    ax.add_geometries([track], ccrs.PlateCarree(),
                      facecolor='none',edgecolor=edgecolor,linewidth=4)


    # Plot a buffer around the track
    if buff:
        track_buffer=track.buffer(buff[0])

        ax.add_geometries([track_buffer], ccrs.PlateCarree(),
                          facecolor='#C8A2C8', alpha=0.5)


    if match_date_to:
        return xax,pmins[ifc_diff:ifc_last+1]



def tc_track(data_path):
    "Open an ascii file containing TC track and intensity information"

    with open(data_path,'r') as fp:

        # Skip header
        next(fp)

        # Make an empty list
        lats=[]
        lons=[]

        for line in fp:
            # Get lat-lon
            lats.append(float(line.split()[4]))
            lons.append(float(line.split()[5]))

        return lats,lons


def tc_depth(data_path):
    "Open an ascii file containing TC track and intensity information"

    with open(data_path,'r') as fp:

        # Skip header
        next(fp)

        # Make an empty list
        dates=[]
        pmins=[]

        for line in fp:
            # Get date elements
            yy=int(line.split()[0])
            mm=int(line.split()[1])
            dd=int(line.split()[2])
            hh=int(line.split()[3])

            dates.append(datetime(yy,mm,dd,hh))

            # Get pressure in hPa, convert Pa
            pmins.append(float(line.split()[6])*100.)

        return dates,pmins



def col_list_data(data_struct,plot_dict,plot_vars=[]):
    "Set colours for tc track plot"

    cols=sns.color_palette(n_colors=len(plot_vars))

    try:
        plot_dict['fig_ens_col']
    except KeyError:
        pass
    else:
        if plot_dict['fig_ens_col']:
            cols=[]
            for ivar in plot_vars:
                # Ensemble member color
                if ivar['ens']=='member':
                    cols.append(plot_dict['fig_ens_col'])

                # Control member color
                elif ivar['ens']=='ctrl':
                    try:
                        plot_dict['fig_ctrl_col']
                    except KeyError:
                        cols.append('k')
                    else:
                        cols.append(plot_dict['fig_ctrl_col'])

                # Ensemble member color
                elif ivar['ens']=='ensmean':
                    try:
                        plot_dict['fig_ensm_col']
                    except KeyError:
                        cols.append('g')
                    else:
                        cols.append(plot_dict['fig_ensm_col'])
            
    return cols
