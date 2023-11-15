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
from typing import Union

import numpy as np
import pandas as pd
import scipy.optimize as optimize

import pytdoa.geodesy as geodesy
from pytdoa.mlat import exact, lls, nlls


################################################################################
# Optimization funtions
def linoptim(
    sensors: Union[pd.DataFrame, np.ndarray],
    tdoas: Union[pd.DataFrame, np.ndarray],
    input_type: str = "xyz",
) -> np.ndarray:
    """
    Obtain the position by linear methods

    Parameters:
    sensors: DataFrame with 'latitude','longitude','height' or 'x','y','z' parameters
    tdoas: array of tdoa values of shape(n-1,1)
    input_type: How positions are represented ([llh | xyz])

    Returns:
    np.array([lat,lon(,altitude)]|[x,y(,z)]) with the resulting coordinates in the desired inputA
    """

    # Sanity checks for input type
    if input_type == "llh":
        reference_c = np.mean(sensors, axis=0)

        sensors_xyz = geodesy.latlon2xy(
            sensors[:, 0], sensors[:, 1], reference_c[0], reference_c[1]
        )
    elif input_type == "xyz":
        sensors_xyz = sensors

    else:
        raise RuntimeError("Unsupported input type. Supported types are: [llh | xyz]")

    # Check whether we call the exact or the least squares routine
    if sensors.shape[0] > 3:
        res = lls.lls(sensors_xyz, tdoas)
    else:
        res = exact.fang(sensors_xyz, tdoas)

    # Output checks (it could happen that we get imaginary or nan results)
    if res.shape[0] == 0:
        print("No solution found")
        res = np.array([[np.inf, np.inf]])
    elif res.shape[0] > 1:
        print("Multiple solutions")
        res = res[0, :]

    if input_type == "llh":
        return geodesy.xy2latlon(
            res[:, 0], res[:, 1], reference_c[0], reference_c[1]
        ).squeeze()
    else:
        return res


def brutefoptim(
    sensors: np.ndarray,
    tdoas: np.ndarray,
    combinations: np.ndarray,
    xrange: float = 200.0,
    yrange: float = 200.0,
    step: float = 5.0,
    epsilon: float = 1e-4,
    maxiter: int = 10,
    workers: int = 1,
    input_type: str = "xyz",
    num_dim: int = 3,
) -> np.ndarray:
    """
    Obtain the position by brute force

    Parameters:
    sensors: Numpy array with 'llh' or 'xyz' parameters
    tdoas: array of tdoa values of shape(n-1,1)
    combinations: list with sensor combinations when computing tdoa pairs
    input_type: How positions are represented ([llh | xyz])
    num_dim: How many dimensions (2D/3D positioning)

    Returns:
    np.array([llh | xyz]) with the resulting latitude and longitude
    """

    if input_type == "llh":
        sensors_xyz = geodesy.llh2ecef(sensors)

    elif input_type == "xyz":
        sensors_xyz = sensors

    else:
        raise RuntimeError("Unsupported input type. Supported types are: [llh | xyz]")

    X0 = np.mean(sensors_xyz, axis=0)

    optimfun = lambda X: nlls.nlls(X, sensors_xyz, tdoas, combinations)

    Xr = np.array([X0[0], X0[1]])
    F_prev = None
    x = xrange
    y = yrange
    st = step
    Ns = int(2 * xrange / step)
    for i in range(maxiter):
        xyrange = (slice(Xr[0] - x, Xr[0] + x, st), slice(Xr[1] - y, Xr[1] + y, st))

        summary = optimize.brute(optimfun, xyrange, full_output=True, workers=workers)

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

    print("Reached maximum number of iterations")
    return Xr


def nonlinoptim(
    sensors: np.ndarray,
    tdoas: np.ndarray,
    combinations: np.ndarray,
    num_dim: int = 3,
    p0: np.ndarray = None,
    input_type: str = "xyz",
    method: str = "L-BFGS-B",
    use_offset: bool = False,
    l: float = .2,
    s: float = 1.0,
) -> np.ndarray:
    """
    Obtain the position by non linear methods

    Parameters:
    sensors: numpy array with 'latitude', 'longitude', 'height' parameters or 'x', 'y', 'z'
    tdoas: array of tdoa values of shape(n-1,1)
    combinations: array with combinations per sensor
    num_dim: Number of dimensions (2/3)
    p0: Initial guess for starting position
    input_type: How positions are represented ([llh | xyz])
    method: Method for minimization
    use_offset: Whether to optimize the offsets of the receivers
    sigma: Standard deviation of the initialization of the offsets

    Returns:
    np.array([llh | xyz]) with the resulting latitude and longitude
    """

    if input_type == "llh":
        sensors_xyz = geodesy.llh2ecef(sensors)
    elif input_type == "xyz":
        sensors_xyz = sensors
    else:
        raise RuntimeError("Unsupported input type. Supported types are: [llh | xyz]")

    sensors_mean = np.mean(sensors_xyz, axis=0)
    sensors_xyz = sensors_xyz - sensors_mean

    # Minimization routine
    # If no initial point is given we start at the center
    if p0 is None:
        X0 = np.zeros(shape=(num_dim, 1))
    else:
        if input_type == "llh":
            X0 = (geodesy.llh2ecef(p0.reshape(-1, num_dim)) - sensors_mean).reshape(
                num_dim, 1
            )
        else:
            X0 = p0 - sensors_mean

    # If user selects to optimize using the offset, we'll have to set a different group of lambdas
    if use_offset:
        rnd_offset = np.zeros(sensors.shape[0])
        X0 = np.append(X0,rnd_offset)
        optimfun = lambda X: nlls.nlls_with_offset(X, sensors_xyz, tdoas, combinations,l=l)
        jac = lambda X: nlls.nlls_with_offset_der(X, sensors_xyz, tdoas, combinations,l=l)
        lb = np.min(sensors_xyz,axis=0) - s
        lb = np.append(lb,-s*np.ones(sensors_xyz.shape[0]))
        ub = np.max(sensors_xyz,axis=0) + s
        ub = np.append(ub,s*np.ones(sensors_xyz.shape[0]))
        bounds = optimize.Bounds(lb, ub, True)
    else:
        optimfun = lambda X: nlls.nlls(X, sensors_xyz, tdoas, combinations)
        jac = lambda X: nlls.nlls_der(X, sensors_xyz, tdoas, combinations)
        
        bounds = optimize.Bounds(np.min(sensors_xyz,axis=0)-s,np.max(sensors_xyz,axis=0)+s, False)

    # Just call the optimization routine now
    summary = optimize.minimize(optimfun, X0, method=method, jac=jac, bounds=bounds)
    res = np.array(summary.x, copy=False)

    if input_type == "llh":
        return geodesy.ecef2llh(res[0:num_dim] + sensors_mean).squeeze()

    elif input_type == "xyz":
        res[0:num_dim] += sensors_mean
        return res
