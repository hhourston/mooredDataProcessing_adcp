# -*- coding: utf-8 -*-

__author__ = 'diwan'

"""
This script contains functions that make two different kinds of plots of ADCP
netCDF-format data. The types of plots are:
    1. North and East current velocities (one plot containing a subplot for each)
    2. Along and cross-shelf current velocities (one plot containing a subplot for each)

"""


import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import math



def resolve_to_alongcross(u_true, v_true, along_angle):
    along_angle = np.deg2rad(along_angle)
    
    u_along = u_true * np.cos(along_angle) + v_true * np.sin(along_angle)
    u_cross = u_true * np.sin(along_angle) - v_true * np.cos(along_angle)
    
    return u_along, u_cross

""" This function is not useful for data that contains outliers
def get_vminvmax(v):
    minmax = math.ceil(np.nanmax(abs(v * 10))) * 0.1 #to get ceiling with 1 decimal place
    return [-minmax, minmax]
"""

def limit_data(data, bad_bins):

    #  Print station name + number
    print(data.attrs['station'], data.attrs['deployment_number'])
    #print(data['LCEWAP01'])

    if data.attrs['orientation'] == 'up':
        print(data.attrs['instrument_depth'])
        bin_depths = float(data.attrs['instrument_depth']) - data.distance
    else:
        bin_depths = data.attrs['instrument_depth'] + data.distance
    print(bin_depths)
    
    
    # data.time should be limited to the data.time with no NA values
    # AS and CS will need to be likewise limited
    # bins must be limited
    # values bin index (the 10) is a random choice; a better method should be used
    for i, x in enumerate(data.LCEWAP01.values[0,10,]):
        if not math.isnan(x):
            new_first = i
            print(new_first)
            break
    
    for j, y in reversed(list(enumerate(data.LCEWAP01.values[0,10,]))):
        if not math.isnan(y):
            new_last = j
            print(new_last)
            break
    
    # Remove bins where surface backscatter occurs
    time_lim = data.time[new_first:new_last]
    
    if bad_bins != 0:
        bin_depths_lim = bin_depths[:-bad_bins]
    
        EW_lim = data.LCEWAP01.values[:, :-bad_bins, new_first:new_last]
        NS_lim = data.LCNSAP01.values[:, :-bad_bins, new_first:new_last]
    else:
        bin_depths_lim = bin_depths
    
        EW_lim = data.LCEWAP01.values[:, :, new_first:new_last]
        NS_lim = data.LCNSAP01.values[:, :, new_first:new_last]
    
    return time_lim, bin_depths_lim, NS_lim, EW_lim



def make_pcolor_ne(data, time_lim, bin_depths_lim, NS_lim, EW_lim):
    print('in ne')
    
    #vminvmax = get_vminvmax(NS_lim)
    vminvmax = [-1.0, 1.0]
    fig = plt.figure(figsize=(13.75, 10))
    ax = fig.add_subplot(2, 1, 1)

    f1 = ax.pcolor(time_lim, bin_depths_lim, NS_lim[0, :, :], cmap='RdBu_r', vmin=vminvmax[0], vmax=vminvmax[1])
    cbar = fig.colorbar(f1, shrink=0.8)
    cbar.set_label('Velocity [m s$^{-1}$]', fontsize=14)
    
    ax.set_ylabel('Depth [m]', fontsize=14)
    ax.set_title('Raw ADCP (North) ' + data.attrs['station'] + '-' + data.attrs['deployment_number'] + ' ' + str(int(data.instrument_depth)) + 'm', fontsize=14)
    if data.orientation == 'up':
        plt.gca().invert_yaxis()

    ax2 = fig.add_subplot(2, 1, 2)
    
    #vminvmax = get_vminvmax(EW_lim)
    
    f2 = ax2.pcolor(time_lim, bin_depths_lim, EW_lim[0, :, :], cmap='RdBu_r', vmin=vminvmax[0], vmax=vminvmax[1])
    cbar = fig.colorbar(f2, shrink=0.8)
    cbar.set_label('Velocity [m s$^{-1}$]', fontsize=14)

    ax2.set_ylabel('Depth [m]', fontsize=14)
    ax2.set_title('Raw ADCP (East) ' + data.attrs['station'] + '-' + data.attrs['deployment_number'] + ' ' + str(int(data.instrument_depth)) + 'm', fontsize=14)
    if data.orientation == 'up':
        plt.gca().invert_yaxis()
    fig.savefig('D:/ADCP_processing/Python_plots/' + data.attrs['station'] + '-' + data.attrs['deployment_number'] + '_{0}m'.format(str(int(data.instrument_depth))) + '-raw_NE_10.png')
    
    return
    
    

def make_pcolor_ac(data, time_lim, bin_depths_lim, NS_lim, EW_lim):
    print('in ac')
    cross_angle = 25 #deg
    along_angle = cross_angle + 90 #deg
    
    CS = np.ones(shape = EW_lim.shape)
    AS = np.ones(shape = NS_lim.shape)
    
    u_along, u_cross = resolve_to_alongcross(EW_lim, NS_lim, along_angle)
    AS = u_along
    CS = u_cross
    
    #vminvmax = get_vminvmax(AS)
    vminvmax = [-1.0, 1.0]
    
    fig = plt.figure(figsize=(13.75, 10))
    
    ax1 = fig.add_subplot(2, 1, 1)
        
    f1 = ax1.pcolor(time_lim, bin_depths_lim, AS[0, :, :], cmap='RdBu_r', vmin=vminvmax[0], vmax=vminvmax[1])
    cbar = fig.colorbar(f1, shrink=0.8)
    cbar.set_label('Velocity [m s$^{-1}$]', fontsize=14)

    ax1.set_ylabel('Depth [m]', fontsize=14)
    ax1.set_title('Raw ADCP (along) ' + str(along_angle) + '$^\circ$ ' + data.attrs['station'] + '-' + data.attrs['deployment_number'] + ' ' + str(math.ceil(data.instrument_depth)) + 'm', fontsize=14)
    if data.orientation == 'up':
        plt.gca().invert_yaxis()
    
    ax2 = fig.add_subplot(2, 1, 2)

    #vminvmax = get_vminvmax(CS)
    
    f2 = ax2.pcolor(time_lim, bin_depths_lim, CS[0, :, :], cmap='RdBu_r', vmin=vminvmax[0], vmax=vminvmax[1])
    cbar = fig.colorbar(f2, shrink=0.8)
    cbar.set_label('Velocity [m s$^{-1}$]', fontsize=14)

    ax2.set_ylabel('Depth [m]', fontsize=14)
    ax2.set_title('Raw ADCP (cross) ' + str(cross_angle) + '$^\circ$ ' + data.attrs['station'] + '-' + data.attrs['deployment_number'] + ' ' + str(math.ceil(data.instrument_depth)) + 'm', fontsize=14)
    if data.orientation == 'up':
        plt.gca().invert_yaxis()
    
    fig.savefig('D:/ADCP_processing/Python_plots/' + data.attrs['station'] + '-' + data.attrs['deployment_number'] + '_{0}m'.format(str(math.ceil(data.instrument_depth))) + '-raw_AC_10.png')
    
    plt.show()
    
    return



#data_dir = 'D:/ADCP_processing/ADCP/'

os.chdir(data_dir)

f = './sample_data/a1_20050503_20050504_0221m.adcp.L1.nc'    
processed_data = xr.open_dataset(f)
bad_bins = 0
    
time_lim, bin_depths_lim, NS_lim, EW_lim = limit_data(processed_data, bad_bins)

make_pcolor_ne(processed_data, time_lim, bin_depths_lim, NS_lim, EW_lim)

make_pcolor_ac(processed_data, time_lim, bin_depths_lim, NS_lim, EW_lim)


#xr.Dataset.close(processed_data)