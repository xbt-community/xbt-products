# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:33:52 2025

@author: Werner
"""

import numpy as np
from scipy.linalg import lstsq

def smooth2d_loess(data, xgrid, ygrid, span_x, span_y, xgrid_est, ygrid_est):
    """
    2D LOESS smoother: smooths a 2D array of data using weighted least squares regression.
    """
    # Check inputs
    if data.ndim != 2:
        raise ValueError("DATA must be a 2D array.")

    ny, nx = data.shape

    if len(xgrid) != nx or len(ygrid) != ny:
        raise ValueError("XGRID or YGRID dimensions do not match DATA dimensions.")

    xgrid = np.asarray(xgrid).flatten()
    ygrid = np.asarray(ygrid).flatten()
    xgrid_est = np.asarray(xgrid_est).flatten()
    ygrid_est = np.asarray(ygrid_est).flatten()

    # Normalize grids by spans
    mxgrid, mygrid = np.meshgrid(xgrid / span_x, ygrid / span_y)
    xgrid_est /= span_x
    ygrid_est /= span_y

    # Preallocate output arrays
    sm_data = np.full((len(ygrid_est), len(xgrid_est)), np.nan)
    flag = np.zeros((len(ygrid_est), len(xgrid_est)), dtype=int)

    # Find locations of NaNs in data
    nan_data_loc = np.isnan(data)

    npt_threshold = 10  # Minimum points required for regression

    for j, ye in enumerate(ygrid_est):
        dy = mygrid - ye
        dy2 = dy ** 2
        jgood = np.where((dy2 < 1) & ~nan_data_loc)

        if len(jgood[0]) >= npt_threshold:
            dy = dy[jgood]
            dy2 = dy2[jgood]
            mxgrid2 = mxgrid[jgood]
            data2 = data[jgood]

            for i, xe in enumerate(xgrid_est):
                dx = mxgrid2 - xe
                dist = dx ** 2 + dy2
                igood = np.where(dist < 1)[0]

                if len(igood) >= npt_threshold:
                    dxsel = dx[igood]
                    dysel = dy[igood]
                    datasel = data2[igood]
                    distsel2 = dist[igood]

                    # Tricubic weighting function
                    w = 1 - distsel2 ** 1.5
                    w = w ** 3

                    # Construct weighted matrix for regression
                    xin = np.column_stack((
                        w, w * dxsel, w * dysel, w * dxsel ** 2, w * dysel ** 2, w * dxsel * dysel
                    ))

                    # Solve weighted least squares
                    B, _, _, _ = lstsq(xin, w * datasel)

                    # Smoothed value
                    sm_value = B[0]
                    sm_data[j, i] = sm_value

                    # Check if out of range
                    if sm_value < datasel.min() or sm_value > datasel.max():
                        flag[j, i] = 1

    return sm_data, flag
