
from __future__ import print_function
import numpy as np
import xarray as xr
import cartopy.crs as ccrs

import configparser
from configparser import SafeConfigParser

from mod_plot import get_varname

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



def data_config(name):

    # Read in configparser configuration file
    parser = SafeConfigParser()
    parser.read("configs/data."+name)

    #for section in parser.sections():        
    #    for name, value in parser.items(section):
    #        print(section, name, value)


    # Return as a dictionary for further use
    my_config_parser_dict = {s:dict(parser.items(s)) for s in parser.sections()}

    # Return sub variables for further use
    main_dict = parse_data_dict(my_config_parser_dict)
    pvars     = parse_vars_dict(my_config_parser_dict)
    operators = parse_oper_dict(my_config_parser_dict)

    return main_dict,pvars,operators



def parse_data_dict(mydict):
    # Form a dictionary similar to what the old main_dict was doing.
    # Parse comma-separated values from the dictionary

    main_dict=[]
    exp=mydict['DATA']
    an=mydict['AN']

    # Read only from DATA part
    for item in exp:
        sub_dict={s:exp[s].split(',') for s in exp}
            
    main_dict.append(sub_dict)


    # Add analysis if defined
    if an:
        for item in an:
            sub_dict={s:an[s].split(',') for s in an}
        
        main_dict.append(sub_dict)

    return main_dict



def parse_vars_dict(mydict):
    # Return variable fields as an array

    vari=[]
    for item in mydict['variables']:
        vari.append(mydict['variables'][item].split(','))
    
    # Change level field into integer number
    for ivar in range(0,len(vari)):
        try:
            vari[ivar][1]
        except IndexError:
            pass
        else:
            vari[ivar][1]=int(vari[ivar][1])

    return vari


def parse_oper_dict(mydict):
    # Return data operators as an array

    oper=[]
    for item in mydict['scores']:
        oper.append(mydict['scores'][item].split(','))

    # Change data source numbers to integers
    for ioper in range(0,len(oper)):
        oper[ioper][1]=int(oper[ioper][1])
        oper[ioper][2]=int(oper[ioper][2])
    
    return oper


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

                    if typ=="ensstd" or typ=="an_test":
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



def create_paths(main_dict):
    "Unroll dictionary elements to construct data paths"

    d_path=[]

    for sub_dict in main_dict:

        # Unroll dictionary elements
        exps=sub_dict['exps']
        dates=sub_dict['dates']
        fnames=sub_dict['types']

        # Construct file paths
        nbpath=0
        for exp in exps:

            basepath=sub_dict['paths'][nbpath]
            if len(sub_dict['paths'])>1:
                nbpath+=1

            for date in dates:
                for fnam in fnames:
                    d_path.append(basepath+exp+"/"+date+"/"+fnam+".nc")

    return d_path



def get_data_layer(pnames,plot_vars,parallel=True):
    "Controls parallel pool or does a serial fetching of data"

    if parallel:
        exit("Not implemented correctly")
        pool=Pool(len(pnames))
        alldata=pool.map(get_data, pnames)
        pool.close()
        pool.join()

    else:
        idata=0
        alldata=[]
        for psource in pnames:
            data=get_data(psource,plot_vars[idata])
            alldata.append(data)

            idata+=1


    return alldata



def get_data(data_path,plot_vars):
    "Open NetCDF file containing data"

    with xr.open_dataset(data_path) as ds:
        ds=xr.open_dataset(data_path)

        # Print out headers
        #print(ds.keys())
        #print 

        # Do a deep copy, can't use the data further down the
        # stream otherwise
        data=xr.Dataset.copy(ds,deep=True)

        # Get variable gribtable name
        item=plot_vars['vars'][0]

        if item=='Z'  : item2='var129'
        if item=='T'  : item2='var130'
        if item=='U'  : item2='var131'
        if item=='V'  : item2='var132'
        if item=='Q'  : item2='var133'
        if item=='VO' : item2='var138'
        if item=='D'  : item2='var155'
            
        if item=='10m gust (3h)'   : item2='var28'
        if item=='10m gust (inst)' : item2='var29'
        if item=='CAPES' : item2='var44'
        if item=='CAPE'  : item2='var59' 
        if item=='TCW'   : item2='var136'
        if item=='MSL'   : item2='var151'
        if item=='TCC'   : item2='var164'
        if item=='U10M'  : item2='var165'
        if item=='V10M'  : item2='var166'
        if item=='T2M'   : item2='var167'
        if item=='D2M'   : item2='var168'
        if item=='TP'    : item2='var228'
        if item=='W10M'  : item2='var255'

        # Pick wanted variable and level (if applicable) from the data
        data_reduced=data[item2]

        if plot_vars['levs']:
            data_reduced=data_reduced.isel(plev=plot_vars['nlevs'][0])

        # Change variable name to ecmwf-gribtable one
        data_reduced.name=item
        
        # Change variable names to more understandable form
        if 1==0:
        #for item in ds.data_vars:
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
        
        return data_reduced




def structure_for_plotting2(data,operators):
    "Unroll data into plottable form"

    idata=0
    data_struct=[]

    # Unroll operators
    for operator in operators:
        
        if operator[0]=='diff':
            data_struct.append(data[operator[2]]-data[operator[1]])

        if operator[0]=='rmse':
            data_struct.append(np.sqrt(np.square(data[operator[2]]-data[operator[1]])))

        if operator[0]=='bias':
            data_struct.append(abs(data[operator[2]]-data[operator[1]]))

    return data_struct



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
    loc_dict.update({'minmax' : ""})

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
