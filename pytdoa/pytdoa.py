#!/usr/bin/env python3

#   Copyright (C) IMDEA Networks Institute 2022
#   This program is free software: you can redistribute it and/or modify
#
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see http://www.gnu.org/licenses/.
#
#   Authors: Yago Lizarribar <yago.lizarribar [at] imdea [dot] org>
#


################################################################################
# Imports
import itertools

import numpy as np
import pandas as pd
from scipy.signal import resample
import scipy.optimize as optimize

from pytdoa.geodesy import geodesy
from pytdoa.mlat import exact, lls, nlls
from pytdoa.util import generate_heatmap, generate_hyperbola

from pytdoa.geodesy.geodesy import SPEED_OF_LIGHT


################################################################################
# PYTDOA functions
def linoptim(sensors, tdoas, input_type="llh"):
    """
    Obtain the position by linear methods

    Parameters:
    sensors: DataFrame with 'latitude','longitude','height' or 'x','y','z' parameters
    tdoas: array of tdoa values of shape(n-1,1)
    input_type: How positions are represented ([llh | xyz])

    Returns:
    np.array([lat,lon]) with the resulting latitude and longitude
    """

    if input_type == "llh":
        if isinstance(sensors, pd.DataFrame):
            sensors_llh = sensors[["latitude", "longitude", "height"]].to_numpy()
        else:
            sensors_llh = sensors
        reference_c = np.mean(sensors_llh, axis=0)

        sensors_xyz = geodesy.latlon2xy(
            sensors_llh[:, 0], sensors_llh[:, 1], reference_c[0], reference_c[1]
        )
    elif input_type == "xyz":
        if isinstance(sensors, pd.DataFrame):
            sensors_xyz = sensors[["x", "y"]].to_numpy()
        else:
            sensors_xyz = sensors
    else:
        RuntimeError("Unsupported input type. Supported types are: [llh | xyz]")

    if sensors.shape[0] > 3:
        res = lls.lls(sensors_xyz, tdoas)
    else:
        res = exact.fang(sensors_xyz, tdoas)

    if res.shape[0] > 1:
        logger.warning("Multiple solutions")
        res = res[0,:]

    if input_type == "llh":
        return geodesy.xy2latlon(res[:,0], res[:,1], reference_c[0], reference_c[1]).squeeze()
    else:
        return res


def brutefoptim(
    sensors,
    tdoas,
    combinations,
    xrange=200,
    yrange=200,
    step=5,
    epsilon=1e-4,
    maxiter=10,
    workers=1,
    input_type="llh",
    num_dim=3
):
    """
    Obtain the position by brute force

    Parameters:
    sensors: DataFrame with 'latitude','longitude','height' or 'x','y','z' parameters
    tdoas: array of tdoa values of shape(n-1,1)
    combinations: list with sensor combinations when computing tdoa pairs
    input_type: How positions are represented ([llh | xyz])
    num_dim: How many dimensions (2D/3D positioning)

    Returns:
    np.array([lat,lon, height]) with the resulting latitude and longitude
    """

    if input_type == "llh":
        if isinstance(sensors, pd.DataFrame):
            sensors_xyz = geodesy.llh2ecef(
                sensors[["latitude", "longitude", "height"]].to_numpy()
            )
        else:
            sensors_xyz = geodesy.llh2ecef(
                sensors
            )

    elif input_type == "xyz":
        if isinstance(sensors, pd.DataFrame):
            if num_dim == 3:
                sensors_xyz = sensors[["x", "y", "z"]].to_numpy()
            elif num_dim == 2:
                sensors_xyz = sensors[["x","y"]].to_numpy()
        else:
            sensors_xyz = sensors
    else:
        RuntimeError("Unsupported input type. Supported types are: [llh | xyz]")

    X0 = np.mean(sensors_xyz, axis=0)

    optimfun = lambda X: nlls.nlls(
        X, sensors_xyz, tdoas, combinations
    )

    Xr = np.array([X0[0], X0[1]])
    F_prev = None
    x = xrange
    y = yrange
    st = step
    Ns = int(2 * xrange / step)
    for i in range(maxiter):
        xyrange = (slice(Xr[0] - x, Xr[0] + x, st), slice(Xr[1] - y, Xr[1] + y, st))

        summary = optimize.brute(
            optimfun, xyrange, full_output=True, finish=optimize.minimize, workers=workers
        )

        # We update all the values for the next iteration
        if F_prev is None:
            Xr = summary[0]
            F_prev = summary[1]

            x = x * 0.1
            y = y * 0.1
            st = 2 * x / Ns
        else:
            Xr = summary[0]
            F = summary[1]

            if np.abs((F - F_prev) / F) < epsilon:
                return Xr

            F_prev = F

            x = x * 0.1
            y = y * 0.1
            st = 2 * x / Ns

    logger.warning("Reached maximum number of iterations")
    return Xr


def nonlinoptim(sensors, tdoas, combinations, num_dim=3, p0=None, input_type="llh", method="BFGS"):
    """
    Obtain the position by non linear methods

    Parameters:
    sensors: DataFrame with 'latitude', 'longitude', 'height' parameters or 'x', 'y', 'z'
    tdoas: array of tdoa values of shape(n-1,1)
    combinations: array with combinations per sensor
    num_dim: Number of dimensions (2/3)
    p0: Initial guess for starting position
    input_type: How positions are represented ([llh | xyz])
    method: Method for minimization

    Returns:
    np.array([lat,lon, height]) with the resulting latitude and longitude
    """

    if input_type == "llh":
        if isinstance(sensors, pd.DataFrame):
            sensors_xyz = geodesy.llh2ecef(
                sensors[["latitude", "longitude", "height"]].to_numpy()
            )
        else:
            sensors_xyz = geodesy.llh2ecef(
                sensors
            )

        sensors_mean = np.mean(sensors_xyz, axis=0)
        sensors_xyz = sensors_xyz - sensors_mean
        optimfun = lambda X: nlls.nlls(X, sensors_xyz, tdoas, combinations)

        # Minimization routine
        # If no initial point is given we start at the center
        if p0 is None:
            X0 = np.zeros(shape=(num_dim, 1))
        else:
            X0 = (geodesy.llh2ecef(p0.reshape(-1,num_dim)) - sensors_mean).reshape(num_dim,1)

        jac = lambda X: nlls.nlls_der(X, sensors_xyz, tdoas, combinations)
        summary = optimize.minimize(optimfun, X0, method=method, jac=jac)

        res = np.array(summary.x, copy=False).reshape(-1, 3)
        return geodesy.ecef2llh(res + sensors_mean).squeeze()

    elif input_type == "xyz":
        if isinstance(sensors, pd.DataFrame):
            if num_dim == 3:
                sensors_xyz = sensors[["x","y","z"]].to_numpy()

            elif num_dim == 2:
                sensors_xyz = sensors[["x","y"]].to_numpy()
        else:
            sensors_xyz = sensors

        if p0 is None:
            X0 = np.zeros(shape=(num_dim, 1))
        else:
            X0 = p0

        optimfun = lambda X: nlls.nlls(X, sensors_xyz, tdoas, combinations)
        jac = lambda X: nlls.nlls_der(X, sensors_xyz, tdoas, combinations)
        summary = optimize.minimize(optimfun, X0, method=method, jac=jac)

        res = np.array(summary.x, copy=False).reshape(-1, num_dim)
        return res

    else:
        RuntimeError("Unsupported input type. Supported types are: [llh | xyz]")

