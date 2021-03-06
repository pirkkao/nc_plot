
from __future__ import print_function
import cartopy.crs as ccrs
from datetime import datetime,timedelta

import re

import configparser
from configparser import SafeConfigParser



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



def data_config(name,exptyp,expdir):

    if not exptyp:
        exptyp="data"

    if not expdir:
        expdir="configs"

    # Read in configparser configuration file
    parser = SafeConfigParser()
    parser.read(expdir+"/"+exptyp+"."+name)

    # Return as a dictionary for further use
    my_config_parser_dict = {s:dict(parser.items(s)) for s in parser.sections()}

    # Return sub variables for further use
    main_dict = parse_data_dict(my_config_parser_dict)
    pvars     = parse_vars_dict(my_config_parser_dict)
    operators = parse_oper_dict(my_config_parser_dict,'scores')
    opers     = parse_oper_dict(my_config_parser_dict,'operations')
    savescore = parse_save_dict(my_config_parser_dict)

    # Unroll possible keywords in main_dict
    main_dict = update_main_dict(main_dict)

    # Include (some) plotting options
    plot_dict2 = parse_plot_dict(my_config_parser_dict)

    
    savescore['fnames']=save_score_names(main_dict,operators,pvars,savescore)

    print("\nMAIN_DICT:"   ,"\n",main_dict)
    print("\nPVARS:"       ,"\n",pvars)
    print("\nOPERATORS:"   ,"\n",opers)
    print("\nSCORES:"      ,"\n",operators)
    print("\nSCORE SAVING:","\n",savescore)
    print("\nPLOTTING:"    ,"\n",plot_dict2,"\n")

    return main_dict,pvars,operators,opers,savescore,plot_dict2



def save_score_names(main_dict,operators,pvars,savescore):
    "Construct file names for both saving and opening operations"

    fname=[]
    if savescore['save_scores'] or savescore['load_scores']:
        
        #for sub_dict in main_dict:
        if 1==1:
            sub_dict=main_dict[0]

            try:
                sub_dict['expand_names']
            except KeyError:
                # AUTO COMPLETE SNAMES
                # Loop over experiment short names
                for exp in sub_dict['snames']:

                    if sub_dict['types'][0]=='ctrl':
                        N='rmse_ctrl'

                        # Loop over dates
                        for date in sub_dict['dates']:

                            # Loop over variables
                            for var in pvars:
                                vtype=var[0]
                                lev=var[1]

                                fname.append(exp+"_"+date+"_"+N+\
                                             "_"+vtype+str(lev)+".nc")

                    else:
                        # Loop over skill score operators
                        ioper=0
                        for operator in operators:
                            otype=operator[0]

                            try:
                                len(operator[2])

                            except TypeError:
                                # Names for RMSE and spread files
                                ftype=sub_dict['types'][ioper]
                                if ftype=='p000':
                                    N=otype+'_ctrl'
                                elif ftype=='ensstd':
                                    N=otype
                                elif ftype=='ensmean':
                                    N='rmse_ensmean'
                                else:
                                    N=otype+"_"+ftype

                                ioper+=1

                            else:
                                N=otype+"_N"+str(len(operator[2]))


                            # Loop over dates
                            for date in sub_dict['dates']:

                                # Loop over variables
                                for var in pvars:
                                    vtype=var[0]
                                    lev=var[1]

                                    fname.append(exp+"_"+date+"_"+N+\
                                                 "_"+vtype+str(lev)+".nc")

            else:
                # EXPAND SNAMES "MANUALLY"
                itype=0
                for exp in sub_dict['snames']:

                    N=sub_dict['types'][itype]
                    itype+=1

                    # Loop over dates
                    for date in sub_dict['dates']:

                        # Loop over variables
                        for var in pvars:
                            vtype=var[0]
                            lev=var[1]

                            fname.append(exp+"_"+date+"_"+N+\
                                                 "_"+vtype+str(lev)+".nc")
                    

    try:
        savescore['fnames']
    except KeyError:
        pass
    else:
        if savescore['fnames']!='default':
            fname=savescore['fnames']

    return fname



def parse_plot_dict(mydict):

    # Parse general plotting 
    #
    tmp=mydict['plot']
    sub_dict={s:tmp[s] for s in tmp}

    # Eval lonlat if not "global"
    try:
        sub_dict['lonlat']!="global"
    except KeyError:
        pass
    else:
        if sub_dict['lonlat']!="global":
            sub_dict['lonlat']=eval(sub_dict['lonlat'])


    # Special treatment for forecast length
    try: 
        tmp['fcsteps']
    except KeyError:
        pass
    else:
        try:
            tmp['fcsteps'].split('/')[1]
        except IndexError:
            sub_dict['fcsteps']=eval(tmp['fcsteps'])
        else:
            tmp_step=tmp['fcsteps'].split('/')

            isteps=[]

            step=int(tmp_step[0])
            last_step=int(tmp_step[2])

            while step <= last_step:

                # Construct indexes
                isteps.append(int(step/eval(tmp['data_fcstep_len'])))

                step=step+int(tmp_step[4])

            sub_dict['fcsteps']=isteps


    # Figure size, nrows, ncols, contour levels
    for key in ['fig_size','fig_nrow','fig_ncol','fig_cf_levs','fig_c_levs','fig_markers',\
                'fig_ens_show','fig_ens_buff']:
        try:
            tmp[key]
        except KeyError:
            pass
        else:
            try:
                eval(tmp[key])
            except SyntaxError:
                sub_dict[key]=False
            else:
                sub_dict[key]=eval(tmp[key])


    # Parse score_plot
    #
    score_plot = mydict['score_plot']

    sub_dict.update({s:score_plot[s].split(',') for s in score_plot})


    # Special treatment for time range
    try:
        score_plot['time'].split(';')[1]
    except KeyError:
        pass
    else:
        time=[]
        time_ranges=score_plot['time'].split(';')
        for itime in range(0,len(time_ranges)):
            time.append(eval(time_ranges[itime]))
        
        sub_dict['time']=time


    # Special treatment for CRPS plotting and time mean
    for key in ["crps","crps_detailed","time_mean","time_mean_load"]:
        try:
            score_plot[key]
        except KeyError:
            sub_dict[key]=False
        else:
            sub_dict[key]=eval(score_plot[key])

    return sub_dict



def parse_save_dict(mydict):
    "Return switches as python logicals"
    
    conf=mydict['save_data']

    savescore={}

    for s in conf:
        try:
            eval(conf[s])
        except NameError:
            savescore[s]=conf[s]
        except SyntaxError:
            savescore[s]=conf[s]
        else:
            savescore[s]=eval(conf[s])

    try:
        savescore['load_scores']
    except KeyError:
        savescore['load_scores']=False

    try:
        savescore['save_scores']
    except KeyError:
        savescore['save_scores']=False
        
    if savescore['load_scores']:
        savescore['save_scores']=False


    return savescore



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

    try:
        sub_dict['dates'][0].split('/')[1]
    except IndexError:
        pass
    else:

        dates=sub_dict['dates'][0].split('/')

        dd  = datetime.strptime(dates[0], "%Y%m%d%H")
        dt1 = datetime.strptime(dates[2], "%Y%m%d%H")

        ddates=[]
        while dd <= dt1:
            ddates.append(dd.strftime("%Y%m%d%H"))
            dd = dd + timedelta(hours=eval(dates[4])*24)

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



def parse_oper_dict(mydict,otype):
    # Return data operators as an array

    oper=[]
    for item in mydict[otype]:
        oper.append(mydict[otype][item].split(','))

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
    
        # Transmute hours to fcsteps
        if oper[ioper][0]=='time':
            fclen=int(mydict['plot']['data_fcstep_len'])
            oper[ioper][2]=int(oper[ioper][2]/fclen)

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

                    # Construct a wildcard for ensmean and std
                    regex_mean=re.compile("ensmea.")
                    regex_std =re.compile("ensst.")

                    if re.match(regex_std,typ) or typ=="an_test":
                        typ_ens.append("ensstd")
                    elif re.match(regex_mean,typ):
                        typ_ens.append("ensmean")
                    elif typ=="ctrl" or typ=="p000":
                        typ_ens.append("ctrl")
                        # MOD FOR WIND BARBS
                        #typ_ens.append("ensmean")
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



def create_paths(main_dict,plot_vars):
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
                    # Iterate over variables if a single input source open
                    # UNLESS AN or CTRL
                    if len(fnames)==1 and len(dates)==1 and fnam!="an_pl" and fnam!="ctrl":
                    # MOD FOR WIND BARBS!
                    #if len(fnames)==1 and len(dates)==1 and fnam!="an_pl":
                        for pvar in plot_vars:
                            d_path.append(basepath+exp+"/"+date+"/"+fnam+".nc")
                            
                    # Else just iterate over fnames and dates
                    else:
                        d_path.append(basepath+exp+"/"+date+"/"+fnam+".nc")

    #print(d_path)
    return d_path



def create_plot_dict():

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
        'fig_nrow'   : [],
        'fig_ncol'   : [],

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

    # Get whole forecast len if requested
    if plot_dict['fcsteps']=="all":
        plot_dict['fcsteps']=range(0,dd[0].sizes['time'])

    return plot_dict



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

    loc_dict.update({'fig_markers'  : False})

    loc_dict.update({'fig_ens_predef': False})
    loc_dict.update({'fig_ens_show'  : False})
    loc_dict.update({'fig_ens_col'   : []})
    loc_dict.update({'fig_ens_buff'   : []})
    loc_dict.update({'fig_ens_alpha'  : []})

    loc_dict.update({'fig_obs_track': False})
    loc_dict.update({'fig_obs_buff':[]})
    loc_dict.update({'fig_obs_match_time': False})

    loc_dict.update({'fig_features': False})
    loc_dict.update({'fig_plot_all_minima': False})

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
                         
                     'fig_markers',\

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
                         
                     'fig_features',\
                     'fig_plot_all_minima']:


        # Replace defaults if a value is given in plot_dict
        try:
            plot_dict[item]
        except KeyError:
            loc_dict.update({item:[]})
        else:
            if plot_dict[item]:
                loc_dict.update({item:plot_dict[item]})

    return loc_dict

