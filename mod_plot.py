
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

# import from mod_tctrack (change to util?)
from mod_tctrack import tc_plot

# import from local
from mod_util import plot_legend, plot_features


# Don't display FutureWarnings
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")



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



def plot_scores(time,data_struct,plot_dict,plot_vars,minmax):
    "Plot scores values as a function of forecast lead time"


    # Create a figure
    fig,ax=plt.subplots(nrows=1,squeeze=0,figsize=plot_dict['fig_size'])

    # Fix the axis handle to be simply ax[0]
    ax=fix_ax(ax)

    plt.tight_layout()
    
    # Loop over data
    for data in data_struct:
        call_score(ax[0],data)

    plot_legend(["b","orange","g"],["CTRL","ENSM","SPREAD"],bbox_loc=(0.3,1.,0,0))


def call_score(ax,data):
    
    dd=data.mean(['lon','lat'])

    print(dd)

    #ax.plot(dd)
    #dd.plot()
    xr.plot.line(dd,ax=ax)



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

            # TEMP SOLUTION FOR SCORES
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


