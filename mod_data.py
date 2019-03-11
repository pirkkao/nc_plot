
from __future__ import print_function
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
from datetime import datetime,timedelta

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
                for imem in range(1,int(sub_dict['nmem'][0])):
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
    "Form a dictionary similar to what the old main_dict was doing.\
     Parse comma-separated values from the dictionary"

    main_dict=[]
    exp=mydict['DATA']
    an=mydict['AN']

    # Read only from DATA part
    sub_dict={s:exp[s].split(',') for s in exp}

    # If dates are given by keywords, generate the wanted dates
    sub_dict=parse_dates_dict(sub_dict)
            
    # Append to main
    main_dict.append(sub_dict)


    # Add analysis if defined
    if an:
        sub_dict={s:an[s].split(',') for s in an}
        
        # If dates are given by keywords, generate the wanted dates
        sub_dict=parse_dates_dict(sub_dict)

        # Append to main
        main_dict.append(sub_dict)


    return main_dict



def parse_dates_dict(sub_dict):
    "Generate date list if keywords 'to' and 'by' are defined"

    if sub_dict['dates'][0].split('/')[1]=='to':

        dates=sub_dict['dates'][0].split('/')

        print("Generate date list:")
        dd  = datetime.strptime(dates[0], "%Y%m%d%H")
        dt1 = datetime.strptime(dates[2], "%Y%m%d%H")

        ddates=[]
        while dd <= dt1:
            ddates.append(dd.strftime("%Y%m%d%H"))
            dd = dd + timedelta(days=int(dates[4]))

        print(ddates)

        sub_dict['dates']=ddates

    return sub_dict



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
        for i in [1,2]:
            try:
                int(oper[ioper][i])
            except ValueError:
                temp=oper[ioper][i].split('/')[0:5:2]
                temp=[int(s) for s in temp]
                oper[ioper][i]=range(temp[0],temp[1],temp[2])
            else:
                oper[ioper][i]=int(oper[ioper][i])
    
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
            for typ in sub_dict['types']:
                for date in sub_dict['dates']:
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
                        
    print()
    print("DATA CONFIGURATION:")
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

            for fnam in fnames:
                for date in dates:
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

        # Do a deep copy, can't use the data further down the
        # stream otherwise [NOT ACTUALLY NEEDED!]
        #data=xr.Dataset.copy(ds,deep=True)

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

        print("GETTING: "+item+" from "+data_path)

        # Pick wanted variable and level (if applicable) from the data
        data_reduced=ds[item2]

        if plot_vars['levs']:
            data_reduced=data_reduced.isel(plev=plot_vars['nlevs'][0])

        # Change variable name to ecmwf-gribtable one
        data_reduced.name=item
        
        
        return data_reduced



def save_score_data(dname,data_struct,nam_list):

    idata=0
    first=True
    for data in data_struct:
        data=data.rename(nam_list[idata])
        if first:
            data.to_netcdf(dname,mode='w')
            first=False
        else:
            data.to_netcdf(dname,mode='a')

        idata+=1



def get_score_data(dname):

    data_struct=[]
    with xr.open_dataset(dname) as ds:
        
        for item in ds:
            data_struct.append(ds[item])
    
    return data_struct



def structure_for_plotting2(data,main_dict,operators):
    "Unroll data into plottable form"

    data_struct=[]

    # Get number of dates
    N = len(main_dict[0]['dates'])
    
    # Get number of data sources
    ndata = len(data)

    # Need to change time axis in order to combine data 
    # from different dates
    for dd in data:
        dd.coords['time'] = range(0,241,6)

    print()

    # Construct name list of operator end results
    nam_list=[]
    nam_ind=0

    # Unroll operators
    for operator in operators:

        print("PROCESSING: "+str(operator))
        
        index=get_index(N,ndata,operator)
        print(" Number of dates",N)
        print(" Number of data sources",ndata)
        print(" Constructed indexes",index)

        # RMSE
        if operator[0]=='rmse':
            data_struct.append(calc_rmse(data,main_dict,index))

        # SPREAD
        if operator[0]=='spread':
            data_struct.append(calc_spread(data,main_dict,index))

        # CRPS
        if operator[0]=='crps':
            crps,fair,crps1,crps2a,crps2b=calc_crps(data,main_dict,index)
            data_struct.append(crps)
            data_struct.append(fair)
            #data_struct.append(crps1)
            #data_struct.append(crps2a)
            #data_struct.append(crps2b)

            nam_list.append("crps"+str(nam_ind))
            nam_list.append("fair"+str(nam_ind))
            nam_ind+=1

    return data_struct, nam_list



def get_index(N,ndata,operator):
    "Check from which element to start reading the data"

    index=[]

    for ii in [1,2]:
        try:
            operator[ii][0]
        except TypeError:

            ipos=operator[ii]
            if ipos == -1:
                ipos = ndata - N
            else:
                ipos*=N

            index.append(ipos)

        else:
            for ipos in operator[ii]:

                if ipos == -1:
                    ipos = ndata - N
                else:
                    ipos*=N

                index.append(ipos)

    return index



def calc_crps(data,main_dict,index):
    "Calculate RMSE of data1 and data2 over N dates"

    # Get number of dates
    N = len(main_dict[0]['dates'])

    # Loop over dates
    crps1=0.
    crps2=0.

    idate=0
    for date in main_dict[0]['dates']:
        
        # i1 should always be AN
        i1=index[0]+idate

        # Loop over ensemble members
        for imem in range(1,len(index)):
            i2=index[imem]+idate

            print(i1,i2,end='\r')
            # Distance to observations
            crps1 = crps1 + np.abs(data[i1]-data[i2])

            # Spread component
            for imem2 in range(1,len(index)):
                i3=index[imem2]+idate
                # Skip calculations with self
                if i3 != i2:
                    print(i1,i2,i3,end='\r')
                    crps2 = crps2 + np.abs(data[i2]-data[i3])

        idate+=1

    M=len(index)-1
    print(N,M)
    crps1 = crps1/M/N
    crps2a = crps2/(2*M**2)/N
    crps2b = crps2/(2*M*(M-1))/N
    

    crps = crps1 - crps2a

    fair = crps1 - crps2b

    return crps,fair,crps1,crps2a,crps2b



def calc_spread(data,main_dict,index):
    "Calculate average spread over N dates"

    # Get number of dates
    N = len(main_dict[0]['dates'])

    # Loop over dates
    spread=0
    idate=0
    for date in main_dict[0]['dates']:
        
        i1=index[0]+idate
        print(i1)

        spread = spread + data[i1]


        idate+=1

    spread = spread/N

    return spread
    


def calc_rmse(data,main_dict,index):
    "Calculate RMSE of data1 and data2 over N dates"

    # Get number of dates
    N = len(main_dict[0]['dates'])

    # Loop over dates
    rmse=0
    idate=0
    for date in main_dict[0]['dates']:
        
        i1=index[0]+idate
        i2=index[1]+idate
        print(i1,i2)

        rmse = rmse + np.square(data[i1]-data[i2])

        #print()
        #print(i1,i2)
        #print(np.square(data[i1]-data[i2]).mean(['lat','lon']))
        #print()

        idate+=1

    rmse = np.sqrt(rmse/N)

    return rmse



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
