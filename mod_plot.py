
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





def plot_data(itime,data_struct,plot_dict,plot_vars,minmax):
    "Plot 2D projections"

    # Generate color maps for data
    cmaps=col_maps_data(data_struct,plot_dict,plot_vars)

    # Create a figure
    fig,ax=create_figure(data_struct,plot_dict)

    # Loop over requested fields
    for idata in range(0,len(data_struct)):
        # Call plotting code layer
        call_plot(ax[idata],data_struct[idata],options=[itime],cmap=cmaps[idata],minmax=minmax[idata])
        
        # Plot additional requested things
        call_plot_add(ax[idata],plot_dict)



def plot_mvar(itime,data_struct,plot_dict,plot_vars,minmax):
    "Plot 2D projections"

    # Generate color maps for data
    cmaps=col_maps_data(data_struct,plot_dict,plot_vars)
    ccont1=[plot_dict['fig_c_col'],plot_dict['fig_c_levs']]
    ccont2=[cmaps[1]              ,plot_dict['fig_cf_levs']]

    # Create a figure
    fig,ax=create_figure(data_struct,plot_dict)

    mim1=minmax[0][0] #min(minmax[0][0],minmax[2][0]) #,minmax[4][0])
    mam1=minmax[0][1] #max(minmax[0][1],minmax[2][1]) #,minmax[4][1])

    mim2=minmax[1][0] #min(minmax[1][0],minmax[3][0]) #,minmax[5][0])
    mam2=minmax[1][1] #max(minmax[1][1],minmax[3][1]) #,minmax[5][1])


    # Call plotting code layer
    call_plot(ax[0],data_struct[0],options=[itime],cmap=ccont1, minmax=[mim1,mam1],plottype='contour')
    call_plot(ax[0],data_struct[1],options=[itime],cmap=ccont2, minmax=[mim2,mam2])
        
    #call_plot(ax[1],data_struct[2],options=[itime],cmap=ccont1, minmax=[mim1,mam1],plottype='contour')
    #call_plot(ax[1],data_struct[3],options=[itime],cmap=ccont2, minmax=[mim2,mam2])

    #call_plot(ax[2],data_struct[4],options=[itime],cmap=ccont1 ,minmax=[mim1,mam1],plottype='contour')
    #call_plot(ax[2],data_struct[5],options=[itime],cmap=ccont2, minmax=[mim2,mam2])

    # Plot additional requested things
    call_plot_add(ax[0],plot_dict)
    #call_plot_add(ax[1],plot_dict)
    #call_plot_add(ax[2],plot_dict)



def call_plot_add(ax,plot_dict):
    "Plot in additional features to the figure"
    
    # Plot Damrey track
    if plot_dict['fig_obs_track']:
        tc_plot(ax,plot_dict['fig_obs_file'],plot_dict['fig_obs_col'],buff=plot_dict['fig_obs_buff'])

    # Plot map features
    if plot_dict['fig_features']:
        plot_features(ax,country=True,lam=True,lonlat=plot_dict['lonlat'])

    # Plot labels if requested
    if plot_dict['fig_title']:
        ax.set_title(plot_dict['fig_title'])
    if plot_dict['fig_ylabel']:
        set_labels(ax,plot_dict['fig_ylabel'])
    if plot_dict['fig_xlabel']:
        set_labels(ax,plot_dict['fig_xlabel'],xory='x')




def create_figure(data_struct,plot_dict):
    "Create figure and return figure and axis handles"

    # Single axis handle
    if (plot_dict['fig_ncol']==1 and plot_dict['fig_nrow']==1):

        # Create a figure
        fig,ax=plt.subplots(nrows=1,squeeze=0,figsize=plot_dict['fig_size'],\
                            subplot_kw={'projection': plot_dict['fig_proj']})


    # Multiple axis handles
    else:
        # Create a figure
        fig,ax=plt.subplots(nrows=plot_dict['fig_nrow'],ncols=plot_dict['fig_ncol'],\
                            figsize=plot_dict['fig_size'],\
                            subplot_kw={'projection': plot_dict['fig_proj']})


    # Fix the axis handle to be simply ax[0]
    ax=fix_ax(ax)

    return fig,ax



def fix_ax(ax):
    "Correct the axis handle when opening a single subplot"

    try:
        ax[0][1]
    except IndexError:
        pass
    except TypeError:
        pass
    else:
        ax_tmp=[]
        for irow in range(0,len(ax)):
            for icol in range(0,len(ax[0])):
                ax_tmp.append(ax[irow][icol])
        return ax_tmp

    try:
        ax[0][0]
    except TypeError:
        pass
    else:
        return ax[0]

    try:
        ax[0]
    except TypeError:
        pass
    else:
        return ax



def call_plot(ax,data,\
              data2=xr.DataArray(data=None),\
              data3=xr.DataArray(data=None),\
              data4=xr.DataArray(data=None),\
              plottype='contourf',options=[],\
              cmap=[],minmax=[],contf=True):
    "Code layer for plotting"

    # Create an empty data structure
    # (multiple input data sources can be used to calculate
    # differences, RMSEs, etc. between the datasets)
    d=[xr.DataArray(data=None),xr.DataArray(data=None),\
       xr.DataArray(data=None),xr.DataArray(data=None)]

    # Choose a timeinstance and level setting from the data
    if get_varname(data)=='TP':
        # Precipitation is cumulative, get the diff of the last fcsteps
        d[0]=data.isel(time=options[0])-data.isel(time=options[0]-1)

    else:
        d[0]=data.isel(time=options[0])

    # Get data information
    #dtime=d[0]['time']

    # Select additional data according to first data set
    if data2.notnull(): d[1]=data2.sel(time=dtime)
    if data3.notnull(): d[2]=data3.sel(time=dtime)
    if data4.notnull(): d[3]=data4.sel(time=dtime)

    #if not minmax.any():
    #    minmax.append(d[0].min().values)
    #    minmax.append(d[0].max().values)

    # Contourf.
    if plottype=="contourf":
        if get_varname(d[0])=='TP':
            contourf_cartopy(ax,d[0],[],[],cmap=cmap)
        else:
            contourf_cartopy(ax,d[0],minmax[0],minmax[1],cmap=cmap)

    if plottype=="contour":
        if get_varname(d[0])=='TP':
            contour_cartopy(ax,d[0],[],[],cmap=cmap)
        else:
            contour_cartopy(ax,d[0],minmax[0],minmax[1],cmap=cmap)
  

def calc_stats(an,data):
    "Calculate global/regional statistics from the given data"

    # GLOBAL NON-WEIGHTED RMSE
    rmse=np.sqrt(((an-data)**2).mean()).values
    #print(np.linalg.norm(an-data)/np.sqrt(160*320))

    return rmse   



def calc_ens_mean(data):
    "Calculate global/regional statistics from the given data"

    ensm=data.mean(dim='time')

    return ensm



def calc_ens_std(data):
    "Calculate global/regional statistics from the given data"

    enss=data.std(dim='time')

    return enss
    

            
def contourf_cartopy(ax,data,fmin,fmax,cmap):
    "Generate a 2D map with cartopy"

    # Cubehelix cmaps are missing len, set it manually
    try:
        len(cmap)
    except TypeError:
        ncolors=10
    else:
        ncolors=len(cmap)-1

    # mvar specific settings for ens spread
    if len(cmap)==2:
        ncolors=cmap[1]
        cmap=cmap[0]

    # Determine contour intervals
    if not fmin==[]:
        conts=np.arange(fmin,fmax,(fmax-fmin)/ncolors)
    else:
        fmin=data.min().values
        fmax=data.max().values
        conts=np.arange(fmin,fmax,(fmax-fmin)/ncolors)

    # Plot
    xr.plot.contourf(data, ax=ax, transform=ccrs.PlateCarree(), \
                     cmap=cmap, levels=conts, extend='max')

            

def contour_cartopy(ax,data,fmin,fmax,cmap):
    "Generate a 2D map with cartopy"

    ccol=cmap[0]
    clen=cmap[1]

    # Determine contour intervals
    if not fmin==[]:
        conts=np.arange(fmin,fmax,(fmax-fmin)/clen)
    else:
        fmin=data.min().values
        fmax=data.max().values
        conts=np.arange(fmin,fmax,(fmax-fmin)/clen)

    # Plot
    cs=xr.plot.contour(data, ax=ax, transform=ccrs.PlateCarree(), \
                     levels=conts,colors=ccol,alpha=0.95)

    ax.clabel(cs,fmt= '%1.0f',fontsize=14)


def plot_features(ax,coast=True,grid=True,country=False,lam=False,lonlat=[]):
    "Plot some cartopy features to the map object"

    # Limited area
    if lonlat=="global":
        ax.set_global()
    elif lonlat:
        ax.set_extent(lonlat)

    # Coast line
    ax.coastlines('50m',edgecolor='gray',facecolor='gray')

    # Grid lines
    gl=ax.gridlines(draw_labels=True, color='gray', alpha=0.4, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = False

    # Draw Vietnamese borders
    plot_borders(ax)



def plot_legend(cols,legend_texts,plot_vars=[],plot_dict=[]):
    "Plot legend box"

    # Default number of legend columns
    ncols=1

    if not legend_texts:
        legend_texts=[]
        for ivar in plot_vars:
            if plot_dict['fig_ens_col']==[] or ivar['ens']!='member':
                legend_texts.append(ivar['legend'])

        bbox_loc=(0.6,1.2,0,0)
        if len(plot_vars) > 18:
            ncols=4
            bbox_loc=(0.9,1.15,0,0)

        elif len(plot_vars) > 12:
            ncols=3
            bbox_loc=(0.8,1.15,0,0)

        elif len(plot_vars) > 6:
            ncols=2
            bbox_loc=(0.7,1.15,0,0)

    else:        
        bbox_loc=[]

    # Create a legend
    legend_artists = [Line([0], [0], color=color, linewidth=2)
                      for color in cols]
    if bbox_loc:
        legend = plt.legend(legend_artists, legend_texts, fancybox=True,
                            loc='upper right', framealpha=0.75,ncol=ncols,
                            bbox_to_anchor=bbox_loc)
    else:
        legend = plt.legend(legend_artists, legend_texts, fancybox=True,
                            loc='upper right', framealpha=0.75,ncol=ncols)

    legend.legendPatch.set_facecolor('wheat')



def plot_borders(ax):
    "Plot Vietnamese borders"

    shpfilename = shapereader.natural_earth(resolution='50m',
                                          category='cultural',
                                          name='admin_0_countries')

    reader = shapereader.Reader(shpfilename)
    countries = reader.records()

    # Find the Vietnam boundary polygon.
    for country in countries:
        if country.attributes['NAME_EN'] == 'Vietnam':
            vietnam = country.geometry
            break
    else:
            raise ValueError('Unable to find the boundary.')

    
    ax.add_geometries([vietnam], ccrs.Geodetic(), edgecolor='black',
                      facecolor='none')




def set_labels(ax,xlabels,xory='y'):
    "Iterate over axis and set labels, need to create a text object for cartopy"

    try:
        len(ax)
    except TypeError:
        ax=[ax]
        xlabels=[xlabels]

    i=0
    for axis in ax:
        if xory=='y':
            axis.text(-0.07, 0.55, xlabels[i], va='bottom', ha='center',
                       rotation='vertical', rotation_mode='anchor',
                       transform=axis.transAxes, fontsize=16)
        else:
            axis.text(0.47, -0.2, xlabels[i], va='bottom', ha='center',
                       rotation='horizontal', rotation_mode='anchor',
                       transform=axis.transAxes, fontsize=16)
        i+=1






def get_varname(data):
    " Find which variable is requested from the data \
     (this is a bit awkward, since xarray is refusing to give \
     the variable name out in taito-setup...)"

    if 'MSL' in str(data):
        vname='MSL'
    elif 'TP' in str(data):
        vname='TP'
    elif 'T2M' in str(data):
        vname='T2M'
    elif 'D2M' in str(data):
        vname='D2M'
    elif 'TCC' in str(data):
        vname='TCC'
    elif 'VO' in str(data):
        vname='VO'
    elif 'U10M' in str(data):
        vname='U10M'
    elif 'V10M' in str(data):
        vname='V10M'
    elif 'W10M' in str(data):
        vname='W10M'
    elif 'U' in str(data):
        vname='U'
    elif 'V' in str(data):
        vname='V'
    elif 'Q' in str(data):
        vname='Q'
    elif 'Z' in str(data):
        vname='Z'
    elif 'D' in str(data):
        vname='D'
    elif 'T' in str(data):
        vname='T'
    else:
        vname='T'

    return vname




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






def col_maps_data(data_struct,plot_dict,plot_vars=[]):
    "Set colormaps for each variable [physical units, standard deviation]"

    cmaps=[]

    clevs=plot_dict['fig_cf_levs']

    idata=0
    for data in data_struct:
        icol=0
        if plot_vars:
            if plot_vars[idata]['ens']=='ensstd':
                icol=1

        cmap=col_maps(get_varname(data),clevs)[icol]

        cmaps.append(cmap)

        idata+=1

    return cmaps




def col_maps(var,clevs):
    "Set colormaps for each variable [physical units, standard deviation]"

    if var=='MSL': 
        cmap=[sns.color_palette("BrBG_r",clevs),sns.cubehelix_palette(n_colors=clevs,start=2.7, light=1, as_cmap=True)]
    elif var=='Z': 
        cmap=[sns.color_palette("BrBG_r",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True, n_colors=clevs)]
    elif var=='T':
        cmap=[sns.color_palette("RdBu_r",clevs),sns.cubehelix_palette(start=2.4, light=1, as_cmap=True, n_colors=clevs)]
    elif var=='Q':
        cmap=[sns.color_palette("PuBu",clevs),sns.cubehelix_palette(start=3.0, light=1, as_cmap=True, n_colors=clevs)]
    elif var=='U' or var=='V':
        cmap=[sns.color_palette("OrRd",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True, n_colors=clevs)]
    elif var=='W10M':
        cmap=[sns.color_palette("cool",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True, n_colors=clevs)]
    elif var=='TP':
        cmap=[sns.color_palette("viridis",clevs-10),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True, n_colors=clevs)]
    else:
        cmap=[sns.color_palette("RdBu_r",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True, n_colors=clevs)]

    return cmap




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
