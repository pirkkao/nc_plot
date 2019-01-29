
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



def update_main_dict(main_dict):
    "Update main_dict for an ensemble"

    tmp_dict=[]
    idata=0
    for sub_dict in main_dict:

        types=[]
        for typ in sub_dict['types']:

            if typ!='ensmemb':
                types.append(typ)
            else:
                for imem in range(1,sub_dict['nmem']):
                    imem=("{:03d}".format(imem))
                    types.append("p"+imem)
        
        tmp_dict.append(sub_dict)
        tmp_dict[idata].update({
            'types' : types,
            })

        idata+=1

    return tmp_dict



def create_vars(pvars,main_dict,plot_dict):
    "Create a dictionary for variables to be plotted"

    ptype=plot_dict['plot_type']
    plegend=plot_dict['fig_legend']

    # Get "length" of main_dict (it can be either appended or
    # contain multiple exps, dates and types).
    # So, unroll all.
    md_sum=0
    md_dates=[]
    md_types=[]
    md_exps=[]
    typ_ens=[]
    for sub_dict in main_dict:
        for exp in sub_dict['exps']:
            for date in sub_dict['dates']:
                for typ in sub_dict['types']:
                    md_sum+=1
                    md_dates.append(date)
                    md_types.append(typ)
                    md_exps.append(exp)

                    if typ=="ensstd":
                        typ_ens.append("ensstd")
                    elif typ=="ensmean":
                        typ_ens.append("ensmean")
                    elif typ=="ctrl" or typ=="p000":
                        typ_ens.append("ctrl")
                    else:
                        typ_ens.append("member")

    if (len(pvars) < md_sum) and len(pvars)>1:
        print("\nNOTE! There's mismatch between variables to be plotted \
and opened data files. Not all data sources will be plotted.\n")

    # Check is plegend a suggestion to form a legend list
    if plegend=='type':
        plegend=[]
        for idata in range(0,md_sum):
            plegend.append(str(md_types[idata]))

    elif plegend=='date':
        plegend=[]
        for idata in range(0,md_sum):
            plegend.append(str(md_dates[idata]))

    elif plegend=='exp':
        plegend=[]
        for idata in range(0,md_sum):
            plegend.append(str(md_exps[idata]))

    elif plegend=='datetype':
        plegend=[]
        for idata in range(0,md_sum):
            plegend.append(str(md_dates[idata])+" "+str(md_types[idata]))

    elif plegend=='expdatetype':
        plegend=[]
        for idata in range(0,md_sum):
            plegend.append(str(md_exps[idata])+" "+str(md_dates[idata])+" "+str(md_types[idata]))


    plot_vars=[]

    # Special treatment for TC track, only MSL needed
    if ptype=="track":
        print("\nOverwriting pvars-list.\n")
        for idata in range(0,md_sum):
            plot_vars.append({
                    'vars'  : ["MSL"],
                    'levs'  : False,
                    'dates' : md_dates[idata],
                    'legend': plegend[idata],
                    'ens'   : typ_ens[idata],
                    })

    else:
        idata=0
        for pvar in pvars:
            # Check is the variable data on vertical levels
            if len(pvar)>1:
                levs=True
                nlevs=pvar[1]
            else:
                levs=False
                nlevs=[]

            plot_vars.append({
                        'vars'  : [pvar[0]],
                        'levs'  : levs,
                        'nlevs' : [nlevs],
                        'dates' : md_dates[idata],
                        'legend': plegend[idata],
                        'ens'   : typ_ens[idata],
                        })
            
            # Only move the iteration forwards if there are more than one 
            # data source open
            if md_sum >1:
                idata+=1

        # If only one variable is given use it for all data sources
        if len(pvars)==1:
            # Start from 1 since plot_vars already contains the 1st dict
            for idata in range(1,md_sum):
                plot_vars.append({
                        'vars'  : [pvar[0]],
                        'levs'  : levs,
                        'nlevs' : [nlevs],
                        'dates' : md_dates[idata],
                        'legend': plegend[idata],
                        'ens'   : typ_ens[idata],
                        })
                        

    print("")
    for plv in plot_vars:
        print(plv)
    print("")

    return plot_vars



def create_paths(main_dict,basepath):
    "Unroll dictionary elements to construct data paths"

    d_path=[]
    for sub_dict in main_dict:

        # Unroll dictionary elements
        exps=sub_dict['exps']
        dates=sub_dict['dates']
        fnames=sub_dict['types']

        # Construct file paths
        for exp in exps:
            for date in dates:
                for fnam in fnames:
                    d_path.append(basepath+exp+"/"+date+"/"+fnam+".nc")

    return d_path



def get_data_layer(pnames,parallel=True):
    "Controls parallel pool or does a serial fetching of data"

    if parallel:
        pool=Pool(len(pnames))
        alldata=pool.map(get_data, pnames)
        pool.close()
        pool.join()

    else:
        alldata=[]
        for psource in pnames:
            alldata.append(get_data(psource))

    return alldata



def get_data(data_path):
    "Open NetCDF file containing data"

    with xr.open_dataset(data_path) as ds:
        ds=xr.open_dataset(data_path)

        # Print out headers
        #print(ds.keys())
        #print 

        # Do a deep copy, can't use the data further down the
        # stream otherwise
        data=xr.Dataset.copy(ds,deep=True)

        # Change variable names to more understandable form
        for item in ds.data_vars:
            if item=='var129': data.rename({'var129': 'Z'},inplace=True)
            if item=='var130': data.rename({'var130': 'T'},inplace=True)
            if item=='var131': data.rename({'var131': 'U'},inplace=True)
            if item=='var132': data.rename({'var132': 'V'},inplace=True)
            if item=='var133': data.rename({'var133': 'Q'},inplace=True)
            if item=='var138': data.rename({'var138': 'VO'},inplace=True)
            if item=='var155': data.rename({'var155': 'D'},inplace=True)
            
            if item=='var28': data.rename({'var28': '10m gust (3h)'},inplace=True)
            if item=='var29': data.rename({'var29': '10m gust (inst)'},inplace=True)
            if item=='var44': data.rename({'var44': 'CAPES'},inplace=True)
            if item=='var59': data.rename({'var59': 'CAPE'},inplace=True)
            if item=='var136': data.rename({'var136': 'TCW'},inplace=True)
            if item=='var151': data.rename({'var151': 'MSL'},inplace=True)
            if item=='var164': data.rename({'var164': 'TCC'},inplace=True)
            if item=='var165': data.rename({'var165': 'U10M'},inplace=True)
            if item=='var166': data.rename({'var166': 'V10M'},inplace=True)
            if item=='var167': data.rename({'var167': 'T2M'},inplace=True)
            if item=='var168': data.rename({'var168': 'D2M'},inplace=True)
            if item=='var228': data.rename({'var228': 'TP'},inplace=True)
            if item=='var255': data.rename({'var255': 'W10M'},inplace=True)
        
        return data


def structure_for_plotting(dd,plot_vars):
    "Unroll data into plottable form"

    idata=0
    data_struct=[]

    # Unroll dictionary and choose data
    for plot_var in plot_vars:
    
        # Loop over variables if requested
        for ivar in plot_var['vars']:
        
            # Loop over levels if requested
            if plot_var['levs']:
                for ilev in plot_var['nlevs']:
                    data_struct.append(dd[idata][ivar].isel(plev=ilev))

            else:
                data_struct.append(dd[idata][ivar])
        
        if len(dd)>1:
            idata+=1

    return data_struct



def configure_plot(plot_dict,plot_vars):
    "Store default plot settings and change them if requested"

    loc_dict={}
    loc_dict.update({'fig_proj':ccrs.PlateCarree()})

    # False and None values are the same, need (probably not)
    # to define False sets here to be included
    loc_dict.update({'fig_title' : False})
    loc_dict.update({'fig_ylabel': False})
    loc_dict.update({'fig_xlabel': False})

    loc_dict.update({'fig_ens_predef': False})
    loc_dict.update({'fig_ens_show'  : False})
    loc_dict.update({'fig_ens_col'   : []})
    loc_dict.update({'fig_ens_buff'   : []})
    loc_dict.update({'fig_ens_alpha'  : []})

    loc_dict.update({'fig_obs_track': False})
    loc_dict.update({'fig_obs_buff':[]})
    loc_dict.update({'fig_obs_match_time': False})

    loc_dict.update({'fig_features': False})

    # Default figure settings for 2dmap
    if plot_dict['plot_type']=='2dmap':
        loc_dict.update({
                'fig_size'    : (14,8),
                'fig_cf_levs' : 40,
                })

    # Default figure settings for mvar
    if plot_dict['plot_type']=='mvar':
        loc_dict.update({
                'fig_size'    : (14,8),
                'fig_nrow'    : 1,
                'fig_ncol'    : 1,
                'fig_cf_levs' : 40,
                'fig_c_levs'  : 30,
                'fig_c_col'   : 'k',
                })

    # Calculate number of columns based on plot_vars length
    if plot_dict['plot_type']!='mvar':
        loc_dict.update({
                'fig_nrow' : len(plot_vars),
                'fig_ncol' : 1,
                })

    # Cycle through all possible fields and add them to dict
    # if not defined
    for item in [    'fcsteps',\
                     'fig_name',\
                     'lonlat',\
                     'minmax',\
                     'plot_type',\
                         
                     'fig_size',\
                     'fig_ncol',\
                     'fig_nrow',\
                         
                     'fig_cf_levs',\
                     'fig_c_levs',\
                     'fig_c_col',\
                         
                     'fig_ens_predef',\
                     'fig_ens_show',\
                     'fig_ens_col',\
                     'fig_ens_buff',\
                     'fig_ens_alpha',\
                     'fig_ctrl_col',\
                     'fig_ensm_col',\
                         
                     'fig_legend',\
                         
                     'fig_title',\
                     'fig_ylabel',\
                     'fig_xlabel',\
                         
                     'fig_proj',\
                         
                     'fig_obs_track',\
                     'fig_obs_file',\
                     'fig_obs_col',\
                     'fig_obs_buff',\
                     'fig_obs_match_time',\
                         
                     'fig_features']:


        # Replace defaults if a value is given in plot_dict
        try:
            plot_dict[item]
        except KeyError:
            loc_dict.update({item:[]})
        else:
            if plot_dict[item]:
                loc_dict.update({item:plot_dict[item]})

    return loc_dict



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



def find_pressure_min(data):
    "Simply find the coordinates of the lowest pressure reading"

    # Find the value
    pmin=data.min().values

    # Filter values with over threshold pressure
    if pmin < 100600.:
        # Find the data point
        ploc=data.where(data==pmin,drop=True)

        return ploc['lon'].values,ploc['lat'].values,pmin

    else:
        return [],[],[]



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
        alpha=0.3
    else:
        alpha=0.9

    if plot_var['ens']!="member" or (plot_var['ens']=="member" and ens_show):
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



def col_maps(var,clevs):
    "Set colormaps for each variable [physical units, standard deviation]"

    if var=='MSL': 
        cmap=[sns.color_palette("BrBG_r",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True)]
    elif var=='Z': 
        cmap=[sns.color_palette("BrBG_r",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True)]
    elif var=='T':
        cmap=[sns.color_palette("RdBu_r",clevs),sns.cubehelix_palette(start=2.4, light=1, as_cmap=True)]
    elif var=='Q':
        cmap=[sns.color_palette("PuBu",clevs),sns.cubehelix_palette(start=3.0, light=1, as_cmap=True)]
    elif var=='U' or var=='V':
        cmap=[sns.color_palette("OrRd",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True)]
    elif var=='W10M':
        cmap=[sns.color_palette("cool",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True)]
    elif var=='TP':
        cmap=[sns.color_palette("viridis",clevs-10),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True)]
    else:
        cmap=[sns.color_palette("RdBu_r",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True)]


    return cmap



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
        if len(plot_vars) > 6:
            ncols=2
            bbox_loc=(0.7,1.2,0,0)

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



def plot_tctracks(data_struct,fcsteps,lonlat):
    "Plot TC tracks for given data"

    # Create a figure
    fig,ax=plt.subplots(nrows=1,squeeze=0,subplot_kw={'projection': ccrs.PlateCarree()})

    # Plot Damrey track
    tc_plot(ax[0][0],'damrey_track.dat','red',buff=[])

    # Plot forecast low track
    cols=['#6a5acd','#00bfff','#0000cd','#00ced1','#5f9ea0','red']

    icol=0
    for data in data_struct:
        pmins=create_tc_track(ax[0][0],data,plot_dict,fcsteps,cols[icol])

        print(pmins)
        icol+=1

    # Legend
    plot_legend(cols,["2017103112","2017110112"])

    # Plot map features
    plot_features(ax[0][0],country=True,lam=True,lonlat=lonlat)



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

    # Plot Damrey track
    if plot_dict['fig_obs_match_time']:
        xax_obs,obs=tc_plot(ax1,'damrey_track.dat','red',buff=plot_dict['fig_obs_buff'],\
                                match_date_to=[dt0,fcsteps])
    else:
        tc_plot(ax1,'damrey_track.dat','red',buff=plot_dict['fig_obs_buff'],\
                    match_date_to=[])

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
    ccont=[plot_dict['fig_c_col'],plot_dict['fig_c_levs']]

    # Create a figure
    fig,ax=create_figure(data_struct,plot_dict)

    mim1=min(minmax[0][0],minmax[2][0],minmax[4][0])
    mam1=max(minmax[0][1],minmax[2][1],minmax[4][1])

    mim2=min(minmax[1][0],minmax[3][0],minmax[5][0])
    mam2=max(minmax[1][1],minmax[3][1],minmax[5][1])


    # Call plotting code layer
    call_plot(ax[0],data_struct[0],options=[itime],cmap=ccont   ,minmax=[mim1,mam1],plottype='contour')
    call_plot(ax[0],data_struct[1],options=[itime],cmap=cmaps[1],minmax=[mim2,mam2])
        
    call_plot(ax[1],data_struct[2],options=[itime],cmap=ccont   ,minmax=[mim1,mam1],plottype='contour')
    call_plot(ax[1],data_struct[3],options=[itime],cmap=cmaps[1],minmax=[mim2,mam2])

    call_plot(ax[2],data_struct[4],options=[itime],cmap=ccont   ,minmax=[mim1,mam1],plottype='contour')
    call_plot(ax[2],data_struct[5],options=[itime],cmap=cmaps[1],minmax=[mim2,mam2])

    # Plot additional requested things
    call_plot_add(ax[0],plot_dict)
    call_plot_add(ax[1],plot_dict)
    call_plot_add(ax[2],plot_dict)
        


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
  
        

def get_minmax_layer(lgetmin,data_struct):
    "Variable level layer for finding min and max of the data, \
    fix min-max to be the same for same variables"

    if lgetmin=="rel":
        minmax=get_minmax(data_struct)

    elif lgetmin=="abs":

        minmax=get_minmax(data_struct)

        # Get variable names to check which elements are the same
        vnames=[]
        for data in data_struct:
            vnames.append(get_varname(data))
        
        # Construct a boolean list telling which name elements are the same
        # I'M SURE THIS COULD'VE BEEN DONE BETTER
        # someone must've written a smart list element comparison function/routine
        bool_list=[]
        vcheck=[]

        for idata in range(0,len(vnames)):
            nf=[]
            for x in vnames:
                if x==vnames[idata]:
                    nf.append(True)
                else:
                    nf.append(False)

            # only save boolean list for the first unique appearance of a variable
            if not vnames[idata] in vcheck:
                vcheck.append(vnames[idata])
                bool_list.append(nf)

        # with the boolean list compare min/max values 
        for ibool in range(0,len(bool_list)):
            # Get the correct elements
            elem=[i for i, x in enumerate(bool_list[ibool]) if x]
            
            # And finally get the abs min
            minmax[elem,0]=min(minmax[elem,0])
            minmax[elem,1]=max(minmax[elem,1])


    elif lgetmin=="":
        minmax=[]
        for data in data_struct:
            minmax.append([[],[]])

    return minmax



def get_minmax(data_struct):
    "Check minimum and maximum of the data"

    mins=[]
    maxs=[]
    for data in data_struct:
        imin=float(data.min().values)
        imax=float(data.max().values)

        # Apply some rounding
        if imin > 1000.:
            imin=round(imin,0)
            imax=round(imax,0)
        elif imin > 100.:
            imin=round(imin,1)
            imax=round(imax,1)
        elif imin >= 0.1:
            imin=round(imin,2)
            imax=round(imax,2)
        elif imin < 0.1: # basically q
            imin=round(imin,5)
            imax=round(imax,5)

        mins.append(imin)
        maxs.append(imax)

    minmax=np.zeros([len(mins),2])
    for i in range(0,len(mins)):
        minmax[i]=[mins[i],maxs[i]]

    return minmax



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
        ncolors=20
    else:
        ncolors=len(cmap)-1

    # Determine contour intervals
    if not fmin==[]:
        conts=np.arange(fmin,fmax,(fmax-fmin)/ncolors)
    else:
        fmin=data.min().values
        fmax=data.max().values
        conts=np.arange(fmin,fmax,(fmax-fmin)/ncolors)

    # Plot
    xr.plot.contourf(data, ax=ax, transform=ccrs.PlateCarree(), \
                     cmap=cmap, levels=conts, extend='both')

            

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
                     levels=conts, extend='both',colors=ccol,alpha=0.4)

    ax.clabel(cs,fmt= '%1.0f',fontsize=14)

