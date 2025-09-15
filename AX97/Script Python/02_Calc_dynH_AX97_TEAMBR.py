#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:18:58 2025

@author: wbarros
"""

import os
import glob
import gsw
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import xarray as xr
import cmocean
import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.io import loadmat, savemat
from scipy.ndimage import uniform_filter1d
from datetime import datetime
import re

########################################################
#               DATA CONFIGURATION
########################################################

# Path to .mat files (output from previous processing)
mat_dir = "/home/wbarros/analisesmovar/MOVAR_analises/netcdfs/arquivos_pro_sal/"

# Create output directory for plots/results
output_dir = "/home/wbarros/analisesmovar/MOVAR_analises/netcdfs/arquivos_pro_sal/plots/"
os.makedirs(output_dir, exist_ok=True)

# Load climatological dataset (Argo)
nc_file = "/data2/tayannepf/toolbox_movar/TS/argo_CLIM_2005-2016.nc"
ref_dep = 500  # Reference depth for climatological comparison
data_str = '202002'  # Reference month/year
month = int(data_str[-2:])

# Read NetCDF climatology
ds = xr.open_dataset(nc_file, decode_times=False)
ds['LONGITUDE'] = (ds['LONGITUDE'] + 360) % 360  # Normalize to 0–360
lon_argo = ds['LONGITUDE'].values
lat_argo = ds['LATITUDE'].values
level = ds['LEVEL'].values
time_argo = ds['TIME'].values + 1
month_index = np.argmin(np.abs(time_argo - month))
level_index = np.argmin(np.abs(level - ref_dep))

# Initialize output lists
sa_list = []
ct_list = []
dyn_height_list = []
all_v_geo = []

# List all .mat files (one per cruise)
mat_files = sorted(glob.glob(os.path.join(mat_dir, "*.mat")))

########################################################
#               PROCESSING LOOP
########################################################

for mat_file in mat_files:
    basename = os.path.splitext(os.path.basename(mat_file))[0]

    # Load temperature/salinity/depth/coords from .mat file
    data = scipy.io.loadmat(mat_file)
    s2 = np.array(data['S2'], dtype=np.float64)
    tt = np.array(data['TT'], dtype=np.float64)
    pp = np.array(data['PP'], dtype=np.float64)
    lon = np.array(data['longitude'], dtype=np.float64).flatten()
    lat = np.array(data['latitude'], dtype=np.float64).flatten()
    p = pp.T  # Pressure in transposed form for GSW functions

    # Reshape coordinates for GSW inputs
    lon2d = lon.reshape(1, -1)
    lat2d = lat.reshape(1, -1)

    # Save original NaN masks
    nan_mask_s2 = np.isnan(s2)
    nan_mask_tt = np.isnan(tt)

    # Compute Absolute Salinity and Conservative Temperature
    sa = gsw.SA_from_SP(s2, p, lon2d, lat2d)
    ct = gsw.CT_from_t(s2, tt, p)
    sigma_theta = gsw.sigma0(sa, ct)

    # Compute dynamic height relative to different reference levels
    dyn_height_100 = gsw.geo_strf_dyn_height(sa, ct, p, 100)
    dyn_height_180 = gsw.geo_strf_dyn_height(sa, ct, p, 180)
    dyn_height_500 = gsw.geo_strf_dyn_height(sa, ct, p, 500)

    # Combine into single dynamic height matrix
    dyn_height_all = np.zeros((81, 49))
    dyn_height_all[:, 0] = dyn_height_100[:, 0]
    dyn_height_all[:, 1] = dyn_height_180[:, 1]
    dyn_height_all[:, 2:] = dyn_height_500[:, 2:]

    ########################################################
    #         MARLOS CORRECTION METHOD (Gpan adjustment)
    ########################################################

    F = dyn_height_all.shape[1]

    # Estimate horizontal spacing (dx) between stations
    if 'lon' in locals():
        dxx = np.abs(np.gradient(lon)) * 111320  # meters
    else:
        dxx = np.ones(F) * 10000  # default 10 km spacing

    gpan = dyn_height_all.copy()

    # Replace zero-values with NaNs
    f = np.where(gpan == 0.0)
    gpan[f] = np.nan

    kp = 50  # Reference level index (zeroed)
    gpan[kp, :] = 0.0

    # Sort gaps by distance from central station
    fi, fj = np.where(np.isnan(gpan))
    lf = len(fi)
    as_ = np.abs(fj - F / 2)
    saa = np.argsort(as_)
    fi = fi[saa]
    fj = fj[saa]

    # Fill missing values using local gradient estimates
    for l in range(lf):
        lj = fj[l]
        lj1 = min(F - 1, lj + 1)
        lj2 = min(F - 1, lj + 2)
        ljm1 = max(0, lj - 1)
        ljm2 = max(0, lj - 2)

        dgpan = gpan[fi[l]:kp + 1, lj2] - gpan[fi[l]:kp + 1, lj1]
        dgpan = 1 / (2 * dxx[lj1]) * dgpan

        signo = -1
        gpan[fi[l]:kp + 1, lj] = gpan[fi[l]:kp + 1, lj1] + signo * dgpan * dxx[lj]

        if np.isnan(gpan[fi[l], lj]):
            signo = 1
            dgpan = gpan[fi[l]:kp + 1, ljm1] - gpan[fi[l]:kp + 1, ljm2]
            dgpan = 1 / (2 * dxx[ljm1]) * dgpan
            gpan[fi[l]:kp + 1, lj] = gpan[fi[l]:kp + 1, ljm1] + signo * dgpan * dxx[lj]

        gpan[:fi[l], lj] += gpan[fi[l], lj]

    ########################################################
    #               ADD CLIMATOLOGICAL OFFSET
    ########################################################

    # Find spatial limits for climatology interpolation
    lon = (lon + 360) % 360
    lonmin_idx = np.argmin(np.abs(lon_argo - lon.min()))
    lonmax_idx = np.argmin(np.abs(lon_argo - lon.max()))
    latmin_idx = np.argmin(np.abs(lat_argo - lat.min()))
    latmax_idx = np.argmin(np.abs(lat_argo - lat.max()))

    if lonmin_idx > lonmax_idx:
        lonmin_idx, lonmax_idx = lonmax_idx, lonmin_idx
    if latmin_idx > latmax_idx:
        latmin_idx, latmax_idx = latmax_idx, latmin_idx

    # Subset and interpolate ADDEP climatology
    addep_sub = ds['ADDEP'].isel(
        LONGITUDE=slice(lonmin_idx, lonmax_idx + 1),
        LATITUDE=slice(latmin_idx, latmax_idx + 1),
        LEVEL=level_index,
        TIME=month_index
    ).values

    addep_sub = np.nan_to_num(addep_sub, nan=np.nan)
    lon_sub = lon_argo[lonmin_idx:lonmax_idx + 1]
    lat_sub = lat_argo[latmin_idx:latmax_idx + 1]
    lon_grid, lat_grid = np.meshgrid(lon_sub, lat_sub)
    lon_flat = lon_grid.ravel()
    lat_flat = lat_grid.ravel()
    addep_flat = addep_sub.T.ravel()
    valid_mask = ~np.isnan(addep_flat) & (addep_flat != 0)

    addep_interp = griddata(
        points=np.array([lon_flat[valid_mask], lat_flat[valid_mask]]).T,
        values=addep_flat[valid_mask],
        xi=np.array([lon, lat]).T,
        method='linear'
    )

    # Fallback to nearest-neighbor if needed
    if np.any(np.isnan(addep_interp)):
        addep_interp_nearest = griddata(
            points=np.array([lon_flat[valid_mask], lat_flat[valid_mask]]).T,
            values=addep_flat[valid_mask],
            xi=np.array([lon, lat]).T,
            method='nearest'
        )
        addep_interp[np.isnan(addep_interp)] = addep_interp_nearest[np.isnan(addep_interp)]

    # Add ADDEP climatology to corrected gpan
    dyn_height_total = gpan + addep_interp[np.newaxis, :]

    # Restore original NaNs
    dyn_height_total[nan_mask_s2 | nan_mask_tt] = np.nan

    ########################################################
    #         FINAL CALCULATIONS: UNITS + VELOCITY
    ########################################################

    dyn_height_total_m = dyn_height_total / 9.81
    dyn_height_total_cm = dyn_height_total_m * 100

    dyn_height_m = dyn_height_all / 9.81
    dyn_height_cm = dyn_height_all / 9.81 * 100

    # Compute geostrophic velocity from corrected dynamic height
    v_geo = gsw.geostrophic_velocity(dyn_height_total, lon, lat)[0]
    v_geo_cm = v_geo * 100

    lon_plot = (lon + 180) % 360 - 180  # Normalize for plotting

    ########################################################
    #               APPEND RESULTS TO LISTS
    ########################################################

    sa_list.append(sa)
    ct_list.append(ct)
    dyn_height_list.append(dyn_height_total)
    all_v_geo.append(v_geo)
