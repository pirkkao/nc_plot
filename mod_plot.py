
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
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from tkinter import *

import re

# import from mod_tctrack (change to util?)
#from mod_tctrack import tc_plot
import mod_tctrack as mtrack

# import from local
from mod_util import plot_legend, plot_features


# Don't display FutureWarnings
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


def plot_master(data_struct,plot_dict,plot_vars,operators,minmax):

    # Open a pdf to plot to
    with PdfPages(plot_dict['fig_name']+'.pdf') as pdf:

        # PLOT 2D MAP
        if plot_dict['plot_type']=="2dmap":

            # Plot each forecast step into its own pdf page
            for fcstep in plot_dict['fcsteps']:

                # Call plotting
                plot_data(fcstep,data_struct,plot_dict,plot_vars,minmax)

                # Save the plot to the pdf and open a new pdf page
                pdf.savefig()
                plt.close()


        # PLOT 2D MAP WITH MULTIPLE VARIABLES
        if plot_dict['plot_type']=="mvar":

            # Plot each forecast step into its own pdf page
            for fcstep in plot_dict['fcsteps']:

                # Call plotting
                plot_mvar(fcstep,data_struct,plot_dict,plot_vars,minmax)

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

                data1=[]
                data2=[]
                
                #if 
                for data in data_struct:

                    print(data.name)

                    # Construct a wildcard for crps and fair
                    regex_crps=re.compile("crps_N.")
                    regex_fair=re.compile("fair_N..")

                    if re.match(regex_crps,data.name):
                        data1.append(data)
                    elif re.match(regex_fair,data.name):
                        data2.append(data)

                    # Or check for exact matches
                    if data.name == 'crps':
                        data1.append(data)
                    elif data.name == 'fair':
                        data2.append(data)
                    elif data.name == 'rmse':
                        data1.append(data)
                    elif data.name == 'spread':
                        data2.append(data)

            elif plot_dict['crps_detailed']:
                data1=data_struct
                data2=[]
            else:
                data1= data_struct[0]
                data2=[]

            nvars=4

            dataX=data1
            dataY=data2

            for ivar in range(0,nvars):
                data1=dataX[ivar:len(dataX)+1:nvars]
                data2=dataY[ivar:len(dataX)+1:nvars]


                print(len(data1),len(data2))

                # Call plotting
                plot_scores3(plot_dict['fcsteps'],data1,plot_dict,plot_vars,minmax)

                # Save the plot to the pdf and open a new pdf page
                pdf.savefig()
                plt.close()

                if data2:
                    # Call plotting
                    plot_scores3(plot_dict['fcsteps'],data2,plot_dict,plot_vars,minmax)

                    # Save the plot to the pdf and open a new pdf page
                    pdf.savefig()
                    plt.close()


        # Plot date by date comparison of CPRS and FAIR

        elif plot_dict['plot_type']=="detail_comparison":
            
            data_crps=[]
            data_fair=[]

            for data in data_struct:
                
                data_crps.append(data[0])
                data_fair.append(data[1])


            plot_date_comparison(data_crps,data_fair,plot_dict,operators,plot_vars)


            # Save the plot to the pdf and open a new pdf page
            pdf.savefig()
            plt.close()



        # PLOT SCORES
        elif plot_dict['plot_type']=="score_rmse":

            data1=[]
            data2=[]
                
            for data in data_struct:

                print(data.name)

                # Or check for exact matches
                if data.name == 'rmse':
                    data1.append(data)
                elif data.name == 'spread':
                    data2.append(data)

            nvars=4

            dataX=data1
            dataY=data2

            for ivar in range(0,nvars):
                data1=dataX[ivar:len(dataX)+1:nvars]
                data2=dataY[ivar:len(dataX)+1:nvars]

                # Call plotting
                plot_scores4(plot_dict['fcsteps'],data1,data2,plot_dict,plot_vars,minmax)


                # Save the plot to the pdf and open a new pdf page
                pdf.savefig()
                plt.close()


        elif plot_dict['plot_type']=="fair_exps":


            data1=[]
                
            for data in data_struct:

                # Construct a wildcard for crps and fair
                if plot_dict['crps']:
                    regex_fair=re.compile("crps_N"+plot_dict['fair_count'][0]+"_..")
                else:
                    regex_fair=re.compile("fair_N"+plot_dict['fair_count'][0]+"_..")

                if re.match(regex_fair,data.name):
                    data1.append(data)


            nvars=4

            dataX=data1

            for ivar in range(0,nvars):
                data1=dataX[ivar:len(dataX)+1:nvars]

                for data in data1:
                    print(data.name)

                plot_scores3(plot_dict['fcsteps'],data1,plot_dict,plot_vars,minmax)

                # Save the plot to the pdf and open a new pdf page
                pdf.savefig()
                plt.close()


        elif plot_dict['plot_type']=="score_diff":
            
            plot_score_comparison_detailed(data_struct,plot_dict,operators,plot_vars)
            
            # Save the plot to the pdf and open a new pdf page
            pdf.savefig()
            plt.close()


        elif plot_dict['plot_type']=="score_diff_avg":
            
            plot_score_comparison(data_struct,plot_dict,operators,plot_vars)
            
            # Save the plot to the pdf and open a new pdf page
            pdf.savefig()
            plt.close()



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
    ccont3=[cmaps[0]              ,plot_dict['fig_cf_levs']]


    # Enlarge fontsizes
    plt.rc('font', size=13)
    plt.rc('axes', labelsize=15)
    plt.rc('xtick',labelsize=15)
    plt.rc('ytick',labelsize=15)

    # Create a figure
    fig,ax=create_figure(data_struct,plot_dict)

    mim1=minmax[0][0] #min(minmax[0][0],minmax[2][0]) #,minmax[4][0])
    mam1=minmax[0][1] #max(minmax[0][1],minmax[2][1]) #,minmax[4][1])

    mim2=minmax[1][0] #min(minmax[1][0],minmax[3][0]) #,minmax[5][0])
    mam2=minmax[1][1] #max(minmax[1][1],minmax[3][1]) #,minmax[5][1])

    # Call plotting code layer
    call_plot(ax[0],data_struct[0],options=[itime],cmap=ccont1, minmax=[mim1,mam1],plottype='contour')
    call_plot(ax[0],data_struct[1],options=[itime],cmap=ccont2, minmax=[mim2,mam2])

    # MOD FOR WIND BARBS
    #call_plot(ax[0],data_struct[0],options=[itime],cmap=ccont3, minmax=[mim1,mam1])
    #call_plot(ax[0],data_struct[1],data2=data_struct[2],options=[itime],cmap=ccont2,minmax=[mim2,mam2],plottype='winds')
        
    #call_plot(ax[1],data_struct[2],options=[itime],cmap=ccont1, minmax=[mim1,mam1],plottype='contour')
    #call_plot(ax[1],data_struct[3],options=[itime],cmap=ccont2, minmax=[mim2,mam2])

    #call_plot(ax[2],data_struct[4],options=[itime],cmap=ccont1 ,minmax=[mim1,mam1],plottype='contour')
    #call_plot(ax[2],data_struct[5],options=[itime],cmap=ccont2, minmax=[mim2,mam2])

    # Plot additional requested things
    call_plot_add(ax[0],plot_dict)
    #call_plot_add(ax[1],plot_dict)
    #call_plot_add(ax[2],plot_dict)

    # Remove border whitespace
    fig.tight_layout()


def plot_scores_crps_vs_fair(time,data_struct,plot_dict,plot_vars,minmax):
    "Plot scores values as a function of forecast lead time"

    print()
    print("CREATING FIGURE")

    # Create a figure
    fig,ax=plt.subplots(nrows=5,ncols=2,figsize=plot_dict['fig_size'])

    # Fix the axis handle to be simply ax[0]
    ax=fix_ax(ax)
    plt.tight_layout()
    
    cols=['r','b']
    styles=['-','--']
    names=['N=8','N=10','N=12','N=20','N=50']

    area='glob'
    time=[range(0,21),range(21,41)]

    for icol in [0,1]:

        for i in range(0,5):
            # Loop over data
            idata=0
            for data in data_struct[i*2:(i+1)*2]:
                call_score(ax[i*2+icol],data,area,time[icol],cols[idata],styles[idata])

                idata+=1


def plot_scores2(time,data_struct,plot_dict,plot_vars,minmax):
    "Plot scores values as a function of forecast lead time"

    print()
    print("CREATING FIGURE")

    areas=['nh','tr','sh']
    time=[range(0,9),range(9,21),range(21,41)]
    time=[range(0,5),range(5,11),range(11,21)]
    ntime=len(time)

    # Create a figure
    fig,ax=plt.subplots(nrows=len(areas),ncols=ntime,figsize=plot_dict['fig_size'])

    # Fix the axis handle to be simply ax[0]
    ax=fix_ax(ax)
    #plt.tight_layout()
    
    cols=['r','b','g','violet','k']
    styles=['-','--',':','-.','-']

    # Loop over FC lead times
    for itime in range(0,ntime):

        # Loop over areas
        for iarea in range(0,len(areas)):

            # Loop over data
            idata=0
            for data in data_struct[0:2]:
                call_score(ax[iarea*ntime+itime],data,areas[iarea],time[itime],cols[0],styles[idata])

                idata+=1

            idata=0
            for data in data_struct[8:10]:
                call_score(ax[iarea*ntime+itime],data,areas[iarea],time[itime],cols[1],styles[idata])

                idata+=1

    #plot_legend(cols[0:5],["8","25","50"],bbox_loc=(0.3,1.,0,0))



def plot_scores3(time,data_struct,plot_dict,plot_vars,minmax):
    "Plot scores values as a function of forecast lead time"

    print()
    print("CREATING FIGURE")

    areas=plot_dict['areas']
    try:
        plot_dict['time']
    except KeyError:
        time=[range(0,41)]
    else:
        time= plot_dict['time']

    ntime=len(time)

    # Enlarge fontsizes
    plt.rc('font', size=13)
    plt.rc('axes', labelsize=15)
    plt.rc('xtick',labelsize=15)
    plt.rc('ytick',labelsize=15)

    # Create a figure
    fig,ax=plt.subplots(nrows=len(areas),ncols=ntime,figsize=plot_dict['fig_size'])

    # Fix the axis handle to be simply ax[0]
    ax=fix_ax(ax)
    #plt.tight_layout()
    
    cols=  plot_dict['cols']
    styles=plot_dict['styles']

    legend_cols= plot_dict['legend_cols']
    legend_names=plot_dict['legend_names']
    legend_styles=plot_dict['legend_styles']

    #for data in data_struct:
    #    print(data.name,data.time)

    # Loop over FC lengths
    for itime in range(0,ntime):

        # Loop over areas
        for iarea in range(0,len(areas)):

            # Loop over data
            idata=0
            for data in data_struct:
                call_score(ax[iarea*ntime+itime],data,areas[iarea],time[itime],cols[idata],styles[idata])

                idata+=1

            # Loop over data
            #idata=0
            #for data in data_struct[1]:
            #    call_score(ax[iarea*ntime+itime],data,areas[iarea],time[itime],cols[idata],styles[idata+5])

            #    idata+=1

            if itime==0 and iarea==0:
                #plot_legend(legend_cols,legend_names,legend_styles,bbox_loc=(-1.45,0.3,0.,0.))
                plot_legend(legend_cols,legend_names,legend_styles,bbox_loc=(0.27,1.,0,0))
                #plot_legend(legend_cols,legend_names,legend_styles,bbox_loc=(0.18,1.,0,0))

            if 1==1:
                ax[iarea*ntime+itime].set_xticks([0,24,48,72,96,120,144,168,192,216,240])
                ax[iarea*ntime+itime].set_xticklabels([0,24,48,72,96,120,144,168,192,216,240])

                #ax[iarea*ntime+itime].set_ylim(0,0.75)

            # LABELS
            ax[iarea*ntime+itime].set_title("")
            if 1==1:
                ax[iarea*ntime+itime].set_xlabel("FORECAST LENGTH IN HOURS")
                ax[iarea*ntime+itime].set_ylabel("RMSE/SPREAD")
                ax[iarea*ntime+itime].set_ylabel("FAIR CPRS")
                #ax[iarea*ntime+itime].set_ylabel("CPRS")
                                

            # YLIMS
            if 1==0:
                ax[iarea*ntime+itime].set_title(areas[iarea])
            #if itime==ntime-1:
            if 1==0:
                if itime == 0:
                    #ax[iarea*ntime+itime].set_ylim(0.2,0.9)
                    #ax[iarea*ntime+itime].set_ylim(0.2,0.95)
                    ax[iarea*ntime+itime].set_ylim(0.2,0.8)
                elif itime == 1:
                    #ax[iarea*ntime+itime].set_ylim(0.85,1.55)
                    #ax[iarea*ntime+itime].set_ylim(0.85,1.7)
                    ax[iarea*ntime+itime].set_ylim(0.65,1.4)
                elif itime == 2:
                    #ax[iarea*ntime+itime].set_ylim(1.5,2.4)
                    #ax[iarea*ntime+itime].set_ylim(1.5,2.55)
                    ax[iarea*ntime+itime].set_ylim(1.15,2.2)

            if 1==0:
                if iarea == 0:
                    ax[iarea*ntime+itime].set_ylim(0.2,0.85)
                elif iarea == 2:
                    ax[iarea*ntime+itime].set_ylim(0.85,1.55)
                elif iarea==1:
                    ax[iarea*ntime+itime].set_ylim(1.55,2.4)
                
    # Remove border whitespace
    fig.tight_layout()



def plot_scores4(time,data_struct,data_struct2,plot_dict,plot_vars,minmax):
    "Plot scores values as a function of forecast lead time"

    print()
    print("CREATING FIGURE")

    areas=plot_dict['areas']
    try:
        plot_dict['time']
    except KeyError:
        time=[range(0,41)]
    else:
        time= plot_dict['time']

    ntime=len(time)
    print(ntime,time,time[0])


    # Enlarge fontsizes
    plt.rc('font', size=13)
    plt.rc('axes', labelsize=15)
    plt.rc('xtick',labelsize=15)
    plt.rc('ytick',labelsize=15)

    # Create a figure
    #fig,ax=plt.subplots(nrows=len(areas),ncols=ntime,figsize=plot_dict['fig_size'])
    fig,ax=plt.subplots(nrows=ntime,ncols=len(areas),figsize=plot_dict['fig_size'])

    # Fix the axis handle to be simply ax[0]
    ax=fix_ax(ax)
    #plt.tight_layout()
    
    cols=  plot_dict['cols']
    styles=plot_dict['styles']

    legend_cols= plot_dict['legend_cols']
    legend_names=plot_dict['legend_names']
    legend_styles=plot_dict['legend_styles']

    #for data in data_struct:
    #    print(data.name,data.time)

    # Loop over FC lengths
    for itime in range(0,ntime):

        # Loop over areas
        for iarea in range(0,len(areas)):

            print(itime,iarea)

            # Loop over data
            idata=0
            for data in data_struct:
                call_score(ax[iarea*ntime+itime],data,areas[iarea],time[itime],cols[idata],styles[idata])

                idata+=1

            # Loop over data
            idata=0
            for data in data_struct2:
                call_score(ax[iarea*ntime+itime],data,areas[iarea],time[itime],cols[idata],'--')

                idata+=1


            # Loop over data
            #idata=0
            #for data in data_struct[1]:
            #    call_score(ax[iarea*ntime+itime],data,areas[iarea],time[itime],cols[idata],styles[idata+5])

            #    idata+=1

            if itime==0 and iarea==0:
                #plot_legend(legend_cols,legend_names,legend_styles,bbox_loc=(0.3,1.,0,0))
                plot_legend(legend_cols,legend_names,legend_styles,bbox_loc=(0.25,1.,0,0))

            ax[iarea*ntime+itime].set_title(areas[iarea])

            if 1==1:
                ax[iarea*ntime+itime].set_xticks([0,24,48,72,96,120,144,168,192,216,240])
                ax[iarea*ntime+itime].set_xticklabels([0,24,48,72,96,120,144,168,192,216,240])


            # LABELS
            ax[iarea*ntime+itime].set_title("")
            if 1==1:
                ax[iarea*ntime+itime].set_xlabel("FORECAST LENGTH IN HOURS")
                ax[iarea*ntime+itime].set_ylabel("RMSE/SPREAD")
                #ax[iarea*ntime+itime].set_ylabel("FAIR CRPS2")

            # YLIMS
            if 1==0:
            #if itime==ntime-1:
                if iarea == 0:
                    ax[iarea*ntime+itime].set_ylim(0.5,3.0)
                elif iarea == 2:
                    ax[iarea*ntime+itime].set_ylim(0.5,4.0)
                elif iarea==1:
                    ax[iarea*ntime+itime].set_ylim(0.2,1.)

    # Remove border whitespace
    fig.tight_layout()



def plot_score_comparison_detailed(data_struct,plot_dict,operators,plot_vars):

    # Create plots
    fig2,ax2=plt.subplots(5,2,figsize=(15,15))

    # Create colors
    cmap=plt.get_cmap('Set1')

    # Setup
    all_varis=[]
    all_times=[]

    # Areas
    areas=[[-90,-20],[-20,20],[20,90]]
    sareas=['SH',     'TR',   'NH']

    # Time ranges
    time_ranges=[[1,8],[8,16],[16,24],[24,32],[32,41]]

    # Construct variable name list
    variable_list=[]
    for var in plot_vars:
        variable_list.append(str(var['vars'][0])+str(var['nlevs'][0]))

    var_size=len(variable_list)

    # Construct ensemble sizes
    sizes=[]
    for oper in operators:
        sizes.append(len(oper[2]))

    # Calculate size of input data 
    data_size=int(len(sizes)*var_size)

    # Calculate number of dates
    date_size=int(len(data_struct)/(data_size))

    print("SIZES")
    print(len(data_struct),len(sizes),var_size,data_size,date_size)


    isize=0
    for M in sizes:

        for idate in range(0,date_size):
            # Set intentation to plots for better distinguishing between variables
            xintent=-0.003

            for ivar in range(0,var_size):
                index = ivar + idate*var_size + isize*var_size*date_size
                print(isize,idate,ivar,index)
                ############
                dd=data_struct[index][0]
                ee=data_struct[index][1]

                #print(data_struct[ivar*len(sizes) + isize][0])

                iarea=0
                for jj in [0,2]:

                    # Do NOT take areal mean yet
                    DD=dd.sel(lat=slice(areas[jj][0],areas[jj][1]))
                    EE=ee.sel(lat=slice(areas[jj][0],areas[jj][1]))

                    # Do division on grid by grid basis
                    FF=DD/EE

                    # Plot against ens size as 1+1/M
                    mval=1+1./M

                    # Setup 1-1-line
                    x=range(0,3)

                    itime=0
                    for time in time_ranges:
                        #GG=FF.mean(['lon','lat'])
                        GG=FF.isel(time=slice(time[0],time[1])).mean(['lon','lat'])

                        HH=GG.mean(['time']).values
                        GG=GG.values

                        xx=[mval+xintent for i in GG.flatten()]
                        XX=[mval         for i in GG.flatten()]

                        xx2=[mval+xintent for i in HH.flatten()]

                        ax2[itime][iarea].plot(x,x,linestyle='-',color='gray',alpha=0.7)
                        #ax2[itime][iarea].scatter(xx,GG.flatten(),color=cmap.colors[ivar],\
                        #                          alpha=0.6)
                        ax2[itime][iarea].scatter(xx2,HH.flatten(),color=cmap.colors[ivar],\
                                                  alpha=0.6)

                        ax2[itime][iarea].set_xlim(1,1.15)
                        ax2[itime][iarea].set_ylim(1,1.24)

                        # Display forecast window
                        ax2[itime][iarea].set_ylabel(str(time[0]*6)+"-"+str(time[1]*6-6),fontsize=16)

                        itime+=1

                    iarea+=1
                xintent+=0.002

        isize+=1



def plot_date_comparison(data_crps,data_fair,plot_dict,operators,plot_vars):


    fig3,ax3=plt.subplots(1,1)

    # FC-WINDOW length
    fclen=41

    # Areas
    areas=[[-90,-20],[-20,20],[20,90]]
    sareas=['SH',     'TR',   'NH']
    
    # Construct variable name list
    variable_list=[]
    for var in plot_vars:
        variable_list.append(str(var['vars'][0])+str(var['nlevs'][0]))

    var_size=len(variable_list)

    # Construct ensemble sizes
    sizes=[]
    for oper in operators:
        sizes.append(len(oper[2]))

    # Calculate size of input data 
    data_size=int(len(sizes)*var_size)

    # Calculate number of dates
    date_size=int(len(data_crps)/(data_size))

    print()
    print("SIZES")
    print(len(data_crps),len(sizes),var_size,data_size,date_size)

    # LOADED DATA ORDER
    # t639_eda+sv_2016120100_crps_N8_T4.nc
    # t639_eda+sv_2016120100_crps_N8_Z3.nc
    # ...
    # t639_eda+sv_2016120900_crps_N8_T4.nc
    # t639_eda+sv_2016120900_crps_N8_Z3.nc
    # ...
    # t639_eda+sv_2016120100_crps_N10_T4.nc
    # t639_eda+sv_2016120100_crps_N10_Z3.nc
    # ...
    # 
    # Size of a single M-set is therefore date_size*var_size

    size_m=date_size*var_size

    index_m=[]
    index_mm=[]

    # CONSTRUCT M for x-axis
    for M in sizes:
        print(M)

        # index_m
        for idate in range(0,date_size):
            xintent=-0.0045
            xintent=-0.009

            for ivar in range(0,var_size):

                # include time
                xxint=-0.002
                for itime in range(0,fclen):
                    index_m.append(1+1./M+xintent+xxint)
                    xxint+=0.0001

                xintent+=0.003
                xintent+=0.006


        # index_mm
        xintent=-0.006
        xintent=-0.009
        index_var=[]

        for ivar in range(0,var_size):

            index_time=[]
            # include time
            xxint=-0.002
            for itime in range(0,fclen):
                index_time.append(1+1./M+xintent+xxint)
                xxint+=0.0001

            index_var.append(index_time)
            xintent+=0.004
            xintent+=0.006

        index_mm.append(index_var)

        
    #print(index_m)

    crps_avg=[]
    fair_avg=[]

    print(data_crps[0])

    for ddata in data_crps:
        #crps_avg.append(ddata.sel(lat=slice(areas[0][0],areas[0][1])).mean(['lon','lat']))
        crps_avg.append(ddata.isel(time=slice(0,42,1)).sel(lat=slice(areas[0][0],areas[0][1])).mean(['lon','lat']))

    for ddata in data_fair:
        #fair_avg.append(ddata.sel(lat=slice(areas[0][0],areas[0][1])).mean(['lon','lat']))
        fair_avg.append(ddata.isel(time=slice(0,42,1)).sel(lat=slice(areas[0][0],areas[0][1])).mean(['lon','lat']))

    print(crps_avg[0])

    #ratio_avg=[]
    #for idata in range(0,len(crps_avg)):
    #    ratio_avg.append(crps_avg[idata]/fair_avg[idata])

    #ax3.scatter(index_m,ratio_avg)

    # DATE MEAN
    crps_avg_dmean=[]
    fair_avg_dmean=[]

    index_date=0
    ratio_fin=[]
    for isize in range(0,len(sizes)):

        ratio_vars=[]
        for ivar in range(0,var_size):
            
            index_dates=[]
            for idate in range(0,date_size):

                # INDEX for dates of same variable and M
                #print(isize,ivar,idate,isize*size_m+ivar+idate*var_size)
                index_dates.append(isize*size_m+ivar+idate*var_size)
                
            ratio_vars.append(np.mean([crps_avg[i]/fair_avg[i] for i in index_dates],axis=0))

        ratio_fin.append(ratio_vars)

    print(index_dates)
    print(ratio_vars)

    #AA=[]
    #for i in index_dates:
    #    AA.append(crps_avg[i])

    #BB=np.mean(AA,axis=0)


    ratio_avg=[]
    for idata in range(0,len(crps_avg)):
        ratio_avg.append(crps_avg[idata]/fair_avg[idata])

    #print(ratio_avg)
    print(len(ratio_avg),len(index_m))

    ax3.scatter(index_m,ratio_avg,color='gray')

    print(len(ratio_fin[0]))
    print(ratio_fin[0])
    print(index_mm[0])

    # Colormap
    colmap=[]
    for i in range(0,var_size):
        colmap.append(cm.get_cmap('PuOr', fclen))

    colmap=cm.get_cmap('PuOr', fclen)

    #print(colmap.shape)
    #print(np.squeeze(colmap).shape)

    for i in range(0,len(sizes)):
        for j in range(0,var_size):
            ax3.scatter(index_mm[i][j],ratio_fin[i][j],color=colmap(range(fclen)),alpha=0.7)


    # Setup 1-1-line
    x=range(0,3)


    ax3.plot(x,x,linestyle='-',color='gray',alpha=0.7)

    #ax3.set_xlim(1,1.15)
    #ax3.set_ylim(1,1.24)
    ax3.set_xlim(1.02,1.08)
    ax3.set_ylim(1,1.15)



def plot_score_comparison(data_struct,plot_dict,operators,plot_vars):


    fig3,ax3=plt.subplots(1,2)

    # Areas
    areas=[[-90,-20],[-20,20],[20,90]]
    sareas=['SH',     'TR',   'NH']

    # Cols
    cols=['c','r','b','g']

    x=range(1,4)
    xmval=[1.,1+1./3.,1+1./2.]
    xmval=[1.,1+1./3.**2,1+1./2.**2]
    xmval=[1.,1+1./3*2.,1+1./2.*2]

    xxx=[1.0,1.5]
    yyy=[1.0,1.35]
    print(x,xmval)

    for data in data_struct:
        print(data.name)

    iarea=0
    for jj in [0,2]:
        ioper=0
        for M in [8,10,12,20,50]:

            # Plot against ens size as 1+1/M
            mval=1+1./M

            allFields=[]

            xint=-0.005

            for ivar in range(0,len(plot_vars)):
            ############
                print(jj,M,ivar,ioper*len(plot_vars),ivar+ioper*(len(plot_vars)+4))
                print(jj,M,ivar,ioper*len(plot_vars),4+ivar+ioper*(len(plot_vars)+4))
                dd=data_struct[0+ivar + ioper*(len(plot_vars)+4)]
                ee=data_struct[4+ivar + ioper*(len(plot_vars)+4)]

                # Average
                #DD=dd.sel(lat=slice(areas[jj][0],areas[jj][1])).mean(['lon','lat'])

                DD=dd.sel(lat=slice(areas[0][0],areas[0][1])).mean(['lon','lat'])

                EE=ee.sel(lat=slice(areas[0][0],areas[0][1])).mean(['lon','lat'])
                #print(ee.coords)
                ##EE=ee.mean(['lon','lat'])

                #EE=EE+ee.sel(lat=slice(areas[2][0],areas[2][1])).mean(['lon','lat'])
                

                # Do division on grid by grid basis
                FF=DD/EE

                #FF=FF.isel(time=slice(0,10))

                allFields.append(FF)

                #xx=[mval for i in allFields[0]]
                xx=[mval+xint for i in FF.values]

                ax3[iarea].scatter(xx,FF.values,color=cols[ivar],\
                                       alpha=0.25)
                   
                xint+=0.0025

            ax3[iarea].plot(x,x,linestyle='-',color='gray',alpha=0.7)
            ax3[iarea].plot(x,xmval,linestyle='--',color='gray',alpha=0.7)
            #ax3[iarea].scatter(xx,allFields[0],color='k',\
            #                       alpha=0.6)


            ioper+=1

        ax3[iarea].set_xlim(1,1.2)
        ax3[iarea].set_ylim(1,1.2)

        # Add titles for some plots
        ax3[iarea].set_title(sareas[jj])

                # Remove ticks from all but bottom
                #if itime < len(ax2)-1:
                #    ax3[iarea].tick_params(labelbottom=False)

                # Remove ticks from all but lefthandside figs
                #if iarea > 0:
                #    ax3[iarea].tick_params(labelleft=False)

        iarea+=1



def call_score(ax,data,area,time,col,style):

    # Cut time if needed
    data=data.isel(time=time)

    # Check are latitudes rolling from North to South or other way
    sign=(data.coords['lat'].values[0]-data.coords['lat'].values[1])

    # Do latitudal slicing of data if requested
    if area=='global':
        True
    elif area=='tr':
        data=data.sel(lat=slice(sign*20,-20*sign))
    elif area=='nh':
        if sign==-1:
            lat1=20
            lat2=80
        else:
            lat1=80
            lat2=20
        data=data.sel(lat=slice(lat1,lat2))
    elif area=='sh':
        if sign==-1:
            lat1=-80
            lat2=-20
        else:
            lat1=-20
            lat2=-80
        data=data.sel(lat=slice(lat1,lat2))

    # Take the areal mean
    dd=data.mean(['lon','lat'])

    xr.plot.line(dd,ax=ax,color=col,linestyle=style,linewidth=2.5)



def call_plot_add(ax,plot_dict):
    "Plot in additional features to the figure"
    
    # Plot Damrey track
    if plot_dict['fig_obs_track']:
        mtrack.tc_plot(ax,plot_dict['fig_obs_file'],plot_dict['fig_obs_col'],buff=plot_dict['fig_obs_buff'],fig_markers=plot_dict['fig_markers'])

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
        ax_tmp=[]
        ax_tmp.append(ax)
        ax=ax_tmp
        return ax
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
    #if get_varname(data)=='TP':
        # Precipitation is cumulative, get the diff of the last fcsteps
    #    d[0]=data.isel(time=options[0])-data.isel(time=options[0]-1)

    #else:
    d[0]=data.isel(time=options[0])

    # Get data information
    #dtime=d[0]['time']

    # Select additional data according to first data set
    if data2.notnull().any(): d[1]=data2.isel(time=options[0])
    if data3.notnull(): d[2]=data3.sel(time=dtime)
    if data4.notnull(): d[3]=data4.sel(time=dtime)


    # Contourf.
    if plottype=="contourf":
        contourf_cartopy(ax,d[0],minmax[0],minmax[1],cmap=cmap)

    if plottype=="contour":
        contour_cartopy(ax, d[0],minmax[0],minmax[1],cmap=cmap)
  
    if plottype=="winds":
        barb_cartopy(ax,d[0],d[1],minmax[0],minmax[1],cmap=cmap)




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

    #conts=[998.,1000.,1002.,1004.,1006.,1008.,1010.,1012.,1014.,1016.,1018.]
    #conts=[0.,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2]
    #conts=[0.,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75]

    # Plot
    xr.plot.contourf(data, ax=ax, transform=ccrs.PlateCarree(), \
                     colors=cmap, levels=conts, extend='both')
                     #colors=cmap, levels=conts, extend='min',cbar_kwargs=dict(label="MSLP"))
                     #colors=cmap, levels=conts, extend='max',cbar_kwargs=dict(label="SDEV(MSLP)"))

            

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
                     levels=conts,colors=ccol,alpha=0.65)

    ax.clabel(cs,fmt= '%1.0f',fontsize=14)

            

def barb_cartopy(ax,data,data2,fmin,fmax,cmap):
    "Generate wind barbs with cartopy"

    ccol=cmap[0]
    clen=cmap[1]

    # Determine contour intervals
    if not fmin==[]:
        conts=np.arange(fmin,fmax,(fmax-fmin)/clen)
    else:
        fmin=data.min().values
        fmax=data.max().values
        conts=np.arange(fmin,fmax,(fmax-fmin)/clen)

    lats=data.lat.data
    lons=data.lon.data

    wind_slice = slice(3, -3, 3)
    # Plot
    ax.barbs(lons[wind_slice],lats[wind_slice],data[wind_slice,wind_slice],data2[wind_slice,wind_slice], \
             transform=ccrs.PlateCarree(),color='b',length=7)

    #ax.clabel(cs,fmt= '%1.0f',fontsize=14)






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

    # Probably working now, remove the above later on

    vname=data.name

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
            #icol=1

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
    elif var=='TESTI':
        cmap=[sns.color_palette("winter",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True, n_colors=clevs)]
    elif var=='MYVARIABLE':
        cmap=[sns.color_palette("winter",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True, n_colors=clevs)]
    else:
        cmap=[sns.color_palette("RdBu_r",clevs),sns.cubehelix_palette(start=2.7, light=1, as_cmap=True, n_colors=clevs)]

    return cmap


