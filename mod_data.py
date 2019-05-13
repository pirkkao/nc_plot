
from __future__ import print_function
import numpy as np
import xarray as xr
import os
import copy

# For parallel excecution
import multiprocessing as mp
import itertools
from functools import partial
from itertools import repeat

# For debugging/performance purposes
import timeit

from mod_plot import get_varname



def get_master(d_path,plot_vars,main_dict,operators,dataoper,savescore,parallel=False):
    "Either open raw nc-files for reading or fetch pre-calculated scores"

    # Open raw variable fields and do skill score calculations if requested
    #
    if not savescore['load_scores']:
        # Fetch all data
        data_struct = get_data_layer(d_path,plot_vars,parallel=parallel)
        print("")

        # Data operations
        if operators:
            data_struct = structure_for_plotting2(data_struct,main_dict,operators,savescore)
            print("")

        if dataoper:
            data_struct = structure_for_plotting3(data_struct,main_dict,dataoper)
            print("")

    # OR load pre-calculated skill score fields
    #
    else:
        data_struct=[]
        for fname in savescore['fnames']:
            data_struct.append(get_score_data(fname,savescore['fpath']))


    # Do temporal averaging if requested
    if savescore['time_mean']:
        data_struct=time_average(data_struct,main_dict,operators)

    return data_struct
    



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
        if item=='TP3'   : item2='var228'
        if item=='TP6'   : item2='var228'
        if item=='TP12'  : item2='var228'
        if item=='W10M'  : item2='var255'

        print("GETTING: "+item+" from "+data_path)

        # Pick wanted variable and level (if applicable) from the data
        data_reduced=ds[item2]

        if plot_vars['levs']:
            # Try whether z-axis is plev or lev
            try:
                data_reduced['plev']
            except KeyError:
                data_reduced=data_reduced.isel(lev=plot_vars['nlevs'][0])
            else:
                data_reduced=data_reduced.isel(plev=plot_vars['nlevs'][0])


        # Change variable name to ecmwf-gribtable one
        data_reduced.name=item


        # Change total accumulated precip to accumated over 3h time window
        if item2=='var228':
            data_reduced=precip_converter(data_reduced,item)
    
        
        return data_reduced



def precip_converter(data,item):
    "Change total accumulated precip to accumated over chosen time window"
    
    #NOTE! The data must be in 3h time steps

    # Copy the field
    dtemp=data

    if item=="TP3":
        fcskip=1
    elif item=="TP6":
        fcskip=2
    elif item=="TP12":
        fcskip=4
    elif item=="TP":
        fcskip=10000

    dd_temp=[]

    first=True
    init=True
    skipper=1
    # Loop over forecast lengths
    for time in dtemp.coords['time'].values:

        # Start when enough hours have accumulated
        if init and skipper < fcskip:
            skipper+=1
            print("SKipping "+str(time))

            if first:
                first=False
                prev_time=time
                prev_time1=time
                prev_time2=time

            continue

        elif init:
            skipper=1
            if first:
                first=False
                prev_time=time
                prev_time1=time
                prev_time2=time

                dd_temp.append(dtemp.sel(time=time))

            init=False
            continue

        
        # Data substraction
        dd_temp.append(data.sel(time=time) - dtemp.sel(time=prev_time).values)

        # Keep book of what is previous step
        if skipper < fcskip:
            prev_time1=time
            prev_time=prev_time2

            skipper+=1

        elif fcskip==1:
            prev_time=time

        else:
            prev_time2=time
            prev_time=prev_time1

            skipper=1

    data=xr.concat([idat for idat in dd_temp],dim='time')

    return data



def save_score_data(dname,fpath,data_struct,nam_list):

    idata=0
    first=True
    for data in data_struct:
        data=data.rename(nam_list[idata])
        if first:
            data.to_netcdf(fpath+dname,mode='w')
            first=False
        else:
            data.to_netcdf(fpath+dname,mode='a')

        idata+=1



def get_score_data(dname,fpath):

    data_struct=[]
    with xr.open_dataset(fpath+dname) as ds:
        
        for item in ds:
            data_struct.append(ds[item])
    
    return data_struct



def structure_for_plotting2(data,main_dict,operators,savescore):
    "Unroll data into plottable form"

    data_struct=[]

    # Get number of dates
    N = len(main_dict[0]['dates'])
    
    # Get number of data sources
    ndata = len(data)

    # Need to change time axis in order to combine data 
    # from different dates
    ntimes=len(data[0].coords['time'])
    idata=0
    for dd in data:
        # Cut AN if FC length shorter than 10d
        if len(dd.coords['time'])>ntimes:
            dd=dd.isel(time=range(0,ntimes))
            
        dd.coords['time'] = range(0,6*ntimes-1,6)

        data[idata]=dd

        idata+=1

    # Initialize name indexing
    nam_ind=0

    # Number of experiments
    nexp=len(main_dict[0]['exps'])

    iexp=0
    # Loop over experiments
    for exp in main_dict[0]['exps']:

        # Unroll operators
        for operator in operators:

            print("\nPROCESSING: "+exp+" "+str(operator))
        
            index=get_index(N,ndata,operator,nexp,iexp)
            print(" Number of dates",N)
            print(" Number of data sources",ndata)

            for date in main_dict[0]['dates']:
                print("  Date is ",date)
                print("  Constructed indexes",index)

                # RMSE
                if operator[0]=='rmse':
                    tmp=[]
                    tmp.append(calc_rmse(data,index))
                    fields=['rmse']

                # SPREAD
                if operator[0]=='spread':
                    tmp=[]
                    tmp.append(calc_spread(data,index))
                    fields=['spread']

                # CRPS
                if operator[0]=='crps':
                    crps,fair,crps1,crps2a,crps2b=calc_crps(data,index)

                    tmp=[crps,fair,crps1,crps2a,crps2b]
                    fields=['crps','fair','crps1','crps2','crsp2b']

                if savescore['fnames']:
                    print("  Save to:",savescore['fpath'],\
                              savescore['fnames'][nam_ind],"\n")

                    save_score_data(savescore['fnames'][nam_ind],savescore['fpath'],\
                                        tmp,fields)

                data_struct.append(tmp)

                # Increment index list
                index=[x+1 for x in index]

                nam_ind+=1

        iexp+=1

    return data_struct



def structure_for_plotting3(data,main_dict,operators):
    "Unroll data into plottable form"

    data_struct=[]

    # Get number of dates
    N = len(main_dict[0]['dates'])
    
    # Get number of data sources
    ndata = len(data)

    # Initialize name indexing
    nam_ind=0

    # Number of experiments
    nexp=len(main_dict[0]['exps'])

    iexp=0
    # Loop over experiments
    for exp in main_dict[0]['exps']:

        # Unroll operators
        for operator in operators:

            print("\nPROCESSING: "+exp+" "+str(operator))
        
            index=get_index(N,ndata,operator,nexp,iexp)
            print(" Number of dates",N)
            print(" Number of data sources",ndata)

            for date in main_dict[0]['dates']:
                print("  Date is ",date)
                print("  Constructed indexes",index)

                # RMSE
                if operator[0]=='rmse':
                    tmp=calc_rmse(data,index)
                    tmp.name=operator[3]

                # SPREAD
                if operator[0]=='spread':
                    tmp=calc_spread(data,index)

                # DIFF
                if operator[0]=='diff':
                    tmp=data[index[0]]-data[index[1]]

                # No operations
                if operator[0]=='none':
                    tmp=data[index[0]]

                # Divide by a constant
                if operator[0]=='div':
                    tmp=data[index[0]]/float(index[1])

                # Substract a constant
                if operator[0]=='sub':
                    tmp=data[index[0]]-index[1]

                data_struct.append(tmp)

                # Increment index list
                index=[x+1 for x in index]

                nam_ind+=1

        iexp+=1

    return data_struct



def get_index(N,ndata,operator,nexp,iexp):
    "Check from which element to start reading the data"

    index=[]

    # Increment experiment index as total opened data sources minus
    # length of analysis fields divided by number of exps
    iexp=int((ndata-N)/nexp)*iexp

    for ii in [1,2]:
        try:
            operator[ii][0]
        except TypeError:

            ipos=operator[ii]
            if ipos == -1:
                ipos = ndata - N
            else:
                ipos=ipos*N + iexp

            index.append(ipos)

        else:
            for ipos in operator[ii]:

                if ipos == -1:
                    ipos = ndata - N
                else:
                    ipos=ipos*N + iexp

                index.append(ipos)

    return index



def calc_crps(data,index):
    "Calculate RMSE of data1 and data2 over N dates"

    crps1=0.
    crps2=0.
        
    # i1 should always be AN!
    i1=index[0]

    # Performance
    start_time = timeit.default_timer()

    if True:
        crps1_time=0.
        crps2_time=0.

        # Loop over ensemble members
        for imem in range(1,len(index)):
            i2=index[imem]

            # Performance
            start_tmp = timeit.default_timer()

            # Distance to observations
            crps1 = crps1 + np.abs(data[i1]-data[i2])

            # Performance
            crps1_time = crps1_time + (timeit.default_timer() - start_tmp)
            start_tmp = timeit.default_timer()

            # Spread component
            for imem2 in range(1+imem,len(index)):
                i3=index[imem2]
                # Skip calculations with self
                if i3 != i2:
                    crps2 = crps2 + 2.*np.abs(data[i2]-data[i3])

            crps2_time = crps2_time + (timeit.default_timer() - start_tmp)

        print( " 1", crps1_time)
        print( " 2", crps2_time)

    # Performance
    end_time1 = timeit.default_timer() - start_time

    # Performance
    start_time = timeit.default_timer()

    if False:
        pool=mp.Pool(processes=ncpus())
        
        # Performance
        start_tmp = timeit.default_timer()

        # Calculate distance to analysis.
        #
        # Construct indexes for ens members
        forecasts = [data[index[x]] for x in range(1,len(index))]

        # Deep copy
        forecasts = pool.map(do_deepcopy,forecasts)

        # Fix data and AN on function call
        anas = [data[i1] for x in range(1,len(index))]

        # Not that much gain from deepcopying these 
        #anas = pool.map(do_deepcopy,anas)

        # Performance
        print(" 0", timeit.default_timer() - start_tmp)
        start_tmp = timeit.default_timer()

        # Call map to iterate over ens members
        crpsA = pool.starmap(calc_crps_distance,zip(anas,forecasts))

        # Performance
        print(" 1", timeit.default_timer() - start_tmp)

        # Release memory
        anas=[]

        # Sum parallel tasks together
        crpsA = sum(crpsA)

        # SPREAD COMPONENT
        # Combinations is listing only unique combinations (x)
        # 0  1  2  3  4  5
        # 1  o  o  o  o  o
        # 2  x  o  o  o  o
        # 3  x  x  o  o  o
        # 4  x  x  x  o  o
        # 5  x  x  x  x  o
        #
        # Since the diagonal values are always 0 and the permuted
        # elements (1,2) vs (2,1) are exactly the same, we just need
        # to multiply each crpsB value by 2 to get to the correct
        # spread score value.

        # Option 1: pool.map
        if True:
            # Performance
            start_tmp = timeit.default_timer()

            crpsC = 0.

            # Calculate distance between ens members
            #
            for imem in range(0,len(forecasts)-1):
                i2=index[imem]
                # Option A, use data-structure. NOTE: loop should be 1,len(index)
                #
                # Fix data on function call
                #func = partial(calc_crps_distance,data[i2])
                # Construct indexes for ens members
                #forecasts1=[data[index[x]] for x in range(1+imem,len(index))]

                # Option B, use forecasts structure
                #
                # Fix data on function call
                func = partial(calc_crps_distance,forecasts[imem])

                # Construct indexes for ens members
                forecasts1 = [forecasts[x] for x in range(imem+1,len(forecasts))]

                crpsB = pool.map(func,forecasts1)  

                # Release memory
                forecasts1=[]

                # Sum parallel tasks together
                crpsB = sum(cb*2. for cb in crpsB)

                crpsC = crpsC + crpsB

            print(" 3", timeit.default_timer() - start_tmp)

        # Option 2: pool.starmap
        if False:
            # Performance
            start_tmp = timeit.default_timer()

            #forecasts=[index[x] for x in range(1,len(index))]
            #forecasts=list(itertools.combinations(forecasts,2))

            forecasts1 = list(itertools.combinations(range(0,len(forecasts)),2))

            forecasts1 = [[forecasts[x],forecasts[y]] for x,y in forecasts1]
            # deepcopying is painfully slow
            #datas = pool.map(do_deepcopy,datas)

            # Call starmap to iterate over ens members
            crpsB = pool.starmap(calc_crps_distance,datas)

            # Release memory
            datas=[]

            # Sum parallel tasks together
            crpsC = sum(cb*2. for cb in crpsB)

            # Performance
            print(" 3", timeit.default_timer() - start_tmp)


    # Performance
    end_time2 = timeit.default_timer() - start_time
    print("\ņ\n PERFORMANCE")
    print(" Time serial",end_time1)
    print(" Time pool",end_time2)

    print("\n Difference between pool and serial")
    #print(" CRPS A",round(sum(sum(sum(crps1.values - crpsA.values))),0),round(sum(sum(sum(crps1.values))),0))
    #print(" CRPS B",round(sum(sum(sum(crps2.values - crpsB.values))),0),round(sum(sum(sum(crps2.values))),0))

    M=len(index)-1
    crps1 = crps1/M
    crps2a = crps2/(2*M**2)
    crps2b = crps2/(2*M*(M-1))

    print()

    crps = crps1 - crps2a

    fair = crps1 - crps2b

    return crps,fair,crps1,crps2a,crps2b



def do_deepcopy(data):

    return copy.deepcopy(data)



def calc_crps_distance(i1,i2):

    #print(i1,i2)
    # Distance to observations
    #return np.abs(data[i1]-data[i2])
    return np.abs(i1-i2)



def calc_spread(data,index):
    "Get correct field from the data_struct"

    i1=index[0]

    #spread =  data[i1]

    return data[i1]
    


def calc_rmse(data,index):
    "Calculate RMSE of data1 and data2 over N dates"

        
    i1=index[0]
    i2=index[1]

    rmse = np.square(data[i1]-data[i2])
    rmse = np.sqrt(rmse)

    return rmse



def ncpus():
    "Check are we running with multiple CPUs and act accordingly"

    try:
        os.environ['SLURM_NTASKS']
    except KeyError:
        parallel=1
    else:
        parallel=int(os.environ['SLURM_NTASKS'])

    print(" NCPUS:",parallel)

    return parallel



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




def get_minmax_layer(plot_dict,data_struct):
    "Variable level layer for finding min and max of the data, \
    fix min-max to be the same for same variables"

    lgetmin=plot_dict['minmax']
    steps=plot_dict['fcsteps']

    if lgetmin=="rel":
        minmax=get_minmax(data_struct,steps)

    elif lgetmin=="abs":

        minmax=get_minmax(data_struct,steps)

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


    elif lgetmin=="nan":
        minmax=[]
        for data in data_struct:
            minmax.append([[],[]])

    return minmax



def get_minmax(data_struct,steps):
    "Check minimum and maximum of the data"

    iimins=[]
    iimaxs=[]
    # Loop over variables
    for data in data_struct:
        mins=[]
        maxs=[]
        # Loop over forecast length indexes
        for step in steps:
            imin=float(data.isel(time=step).min().values)
            imax=float(data.isel(time=step).max().values)

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

        iimins.append(min(mins))
        iimaxs.append(max(maxs))

    minmax=np.zeros([len(data_struct),2])
    for i in range(0,len(data_struct)):
        minmax[i]=[iimins[i],iimaxs[i]]

    return minmax



def time_average(data_struct,main_dict,operators):

    nexp=len(main_dict[0]['exps'])
    ndates=len(main_dict[0]['dates'])
    noper=len(operators)

    dd=[]

    if operators[0][0]=='crps':
        nscores=2
    else:
        nscores=1

    print("\n","CALCULATING TIME MEAN","\n")

    iexp=0
    for exp in main_dict[0]['exps']:

        ioper=0
        for operator in operators:

            # Loop over scores
            for iscore in range(0,nscores):
                print("Averaging...")
                idate=0
                idd=0

                for date in main_dict[0]['dates']:

                    print("INDEX:",ndates*noper*iexp + ndates*ioper + idate,"DATE:",date,"ISCORE:",iscore)
                    idd = idd + data_struct[ndates*noper*iexp + ndates*ioper + idate][iscore]

                    idate+=1

                idd=idd/ndates
            
                dd.append(idd)

            ioper+=1
        iexp+=1

    return dd
