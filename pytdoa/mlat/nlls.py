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

import numpy as np

from pytdoa.geodesy import geodesy


def nlls(X, positions, tdoas, combinations):
    """
    Solve TDOA equations using the Non-Linear Least Squares approach
    The solutions contain the ecef coordinates of the estimated
    transmitter position
    ---
    """

    # Compute all distances to the sensor
    d = geodesy.ecef_distance(positions, X)

    si = combinations[:, 0]
    sj = combinations[:, 1]

    t = (d[si] - d[sj]).reshape(-1, 1)

    err_sq = np.square(t - tdoas)
    F = 0.5 * np.sum(err_sq)

    return F


def nlls_der(X, positions, tdoas, combinations, eps=1e-3):
    """
    Jacobian computation for TDOA problems using Non-Linear Least
    Squares
    ---
    """

    # Compute all distances to the sensor
    n = positions.shape[1]
    d = geodesy.ecef_distance(positions, X)

    si = combinations[:, 0]
    sj = combinations[:, 1]

    t = (d[si] - d[sj]).reshape(-1, 1)

    err = t - tdoas

    J = np.array([0.0, 0.0])
    for i in range(n):
        Jij = err * (
            (X[i] - positions[si, i]) / (d[si] + eps)
            - (X[i] - positions[sj, i]) / (d[sj] + eps)
        ).reshape(-1, 1)
        J[i] = np.sum(Jij)

    return J


def nlls_with_offset(X, positions, tdoas, combinations, l=.1):
    """
    Solve TDOA equations using the Non-Linear Least Squares approach
    The solutions contain the ecef coordinates of the estimated
    transmitter position. This option includes the offset of each sensor.
    ---
    """

    n = positions.shape[1]
    offsets = X[n:] # We assume first offset to be 0
    d = geodesy.ecef_distance(positions, X[0:n])

    si = combinations[:, 0]
    sj = combinations[:, 1]

    t = (d[si] - d[sj]).reshape(-1, 1)
    o = (offsets[si] - offsets[sj]).reshape(-1, 1)

    err_sq = np.square((t + l*o) - tdoas.reshape(-1,1))
    F = 0.5 * np.sum(err_sq)

    return F


def nlls_with_offset_der(X, positions, tdoas, combinations, eps=1e-3, l=.1):
    """
    Solve TDOA equations using the Non-Linear Least Squares approach
    The solutions contain the ecef coordinates of the estimated
    transmitter position. This option includes the offset of each sensor.
    ---
    """

    # Let's get some parameters from the shape of the input
    n = positions.shape[1]
    m = np.max(X.shape) - n
    offsets = X[n:]
    d = geodesy.ecef_distance(positions, X[0:n])

    si = combinations[:, 0]
    sj = combinations[:, 1]

    t = (d[si] - d[sj]).reshape(-1, 1)
    o = (offsets[si] - offsets[sj]).reshape(-1, 1)

    err = (t + l*o) - tdoas.reshape(-1,1)

    J = np.zeros(np.max(X.shape))
    # Jacobian for position derivatives
    for i in range(n):
        Jij = err * (
            (X[i] - positions[si, i]) / (d[si] + eps)
            - (X[i] - positions[sj, i]) / (d[sj] + eps)
        ).reshape(-1, 1)
        J[i] = np.sum(Jij)
    
    # Jacobian for offset derivatives
    for i in range(0,m):
        p = (si == i) # When on the right side, derivative it will be 1*err
        m = (sj == i) # When on the left side, derivative it will be -1*err

        J[i+n] = l*(np.sum(err[p]) - np.sum(err[m]))

    return J
