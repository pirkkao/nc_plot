
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




def plot_legend(cols,legend_texts,plot_vars=[],plot_dict=[],bbox_loc=[]):
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
