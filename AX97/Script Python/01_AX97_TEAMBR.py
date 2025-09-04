#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:16:24 2025

@author: wbarros
"""

import os
import sys
import numpy as np
np.bool = np.bool_  # workaround for deprecated np.bool in recent NumPy
import xarray as xr
from datetime import datetime
import matplotlib.pyplot as plt
import time as pytime
import subprocess
from scipy.io import savemat, loadmat
from scipy.interpolate import griddata
import glob
import scipy.io

# Path to input data (NetCDF files from all cruises)
data_dir = '/home/wbarros/analisesmovar/MOVAR_analises/netcdfs/todos_cruzeiros/'
sys.path.append(data_dir)

# Import custom functions
from smooth2d_loess import smooth2d_loess
from Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe import make_sal

# List all NetCDF files in the directory
nc_files = glob.glob(os.path.join(data_dir, '*.nc'))
n_files = len(nc_files)

# Load longitude, depth, and temperature from each file
data_list = []
for ii in range(n_files):
    ds = xr.open_dataset(nc_files[ii])
    data_list.append({
        'long': ds['longitude'].values,
        'prof': ds['depth'].values,
        'temp': ds['temperature'].values
    })

# Generate standard section grid for AX97 (change if using a different section)
longitude = np.concatenate([
    np.arange(-40.9, -40.0, 0.16667),
    np.arange(-40.0, -29.6, 0.25),
    np.array([-29.65])
])[:, np.newaxis]

# Estimate latitude along the section line
p = np.polyfit([-40.9, -29.65], [-23, -20.55], 1)
latitude = np.polyval(p, longitude)

# Depth levels (0 to 800 m every 10 m)
depth = np.arange(0, 810, 10)
A, B = np.meshgrid(longitude.flatten(), depth)

# Load bathymetric mask for AX97 section
mask_file = "/home/wbarros/analisesmovar/MOVAR_analises/mascara.mat"
mask_data = scipy.io.loadmat(mask_file)
mask1 = mask_data['mask1']

# Loop through all cruises
for i in range(len(data_list)):
    long = data_list[i]['long'].astype(float)
    prof = data_list[i]['prof'].astype(float)
    temp = data_list[i]['temp'].astype(float)

    # Replace 0s with NaN (missing data)
    temp[temp == 0] = np.nan
    prof[prof == 0] = np.nan

    # Fix shape mismatch if needed
    if temp.shape != long.shape:
        print('TEMP and LONG have different shapes, adjusting...')
        long = np.tile(long, (1, temp.shape[1]))

    if temp.shape != prof.shape:
        print('TEMP and DEPTH have different shapes, adjusting...')
        prof = np.tile(prof, (2, temp.shape[0])).T

    # Basic correction for NaNs at edges
    if np.isnan(prof[0, 0]) or np.isnan(prof[0, 3]):
        prof[:, 0] = 0
        prof[0, :] = 0

    # Interpolate temperature onto the standard grid
    method_interp = 'nearest'
    valid = np.where(~np.isnan(temp))
    points = np.column_stack([long[valid], prof[valid]])
    values = temp[valid]
    T12 = griddata(points, values, (A, B), method=method_interp)

    # Smooth the temperature field using LOESS
    T12_float = T12.astype(np.float64)
    lon_flat = longitude.flatten().astype(np.float64)
    depth_float = depth.astype(np.float64)

    T, _ = smooth2d_loess(
        T12_float,
        lon_flat,
        depth_float,
        np.float64(1),  # span in x
        np.float64(50), # span in z
        lon_flat,
        depth_float
    )

    # Apply bathymetry mask
    corrected_mask = np.zeros_like(T)
    corrected_mask[:, :48] = mask1
    T[corrected_mask == 1] = np.nan

    # Prepare variables for salinity estimation
    TT = T.copy()
    PP = np.tile(depth[:, np.newaxis], (1, TT.shape[1]))
    lat = latitude.flatten()
    lon = longitude.flatten()

    # Extract time string from file name (format: YYYYMM)
    time_str = os.path.basename(nc_files[i])[19:-3]
    timet = np.full(len(latitude), float(time_str))

    # Define parameters for salinity model
    Pad = 0
    Pout = depth.copy()
    method_sal = 'Goes'  # Options: 'Goes', 'svd', 'Thacker', etc.

    # Run salinity estimation
    sal = make_sal(TT, PP, lat, lon, timet, pad=Pad, Pout=Pout, method=method_sal)

    print(f'Return type from make_sal: {type(sal)}')

    # Handle cases where make_sal returns a tuple
    if isinstance(sal, tuple):
        print(f'Tuple length: {len(sal)}')
        for idx, item in enumerate(sal):
            print(f'Item {idx} type: {type(item)}')
            if hasattr(item, 'shape'):
                print(f'Item {idx} shape: {item.shape}')
        sal = sal[0]  # Use first element as salinity field

    print(f'Final salinity shape: {sal.shape}')

    # Save results to .mat file
    output_name = f'/home/wbarros/analisesmovar/MOVAR_analises/netcdfs/arquivos_pro_sal/marlos_/{time_str}.mat'
    savemat(output_name, {
        'temperature': TT,
        'sal': sal,
        'time': timet,
        'lat': lat,
        'lon': lon,
        'depth': Pout
    })
