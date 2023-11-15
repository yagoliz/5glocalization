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

from typing import Tuple, Union

import numpy as np
import numpy.typing as npt

# Typing alias for better linting
ndarray_f64 = npt.NDArray[np.float64]
ndarray_i64 = npt.NDArray[np.int64]


def mse(
    x: ndarray_f64,
    y: ndarray_f64,
    tdoas: ndarray_f64,
    rx: ndarray_f64,
    si: ndarray_i64,
    sj: ndarray_i64,
) -> ndarray_f64:
    """
    Evaluate the tdoa cost function at a given x, y vector pair list

    Parameters:
    x: X coordinates as an Nx1-sized array
    y: Y coordinates as an Nx1-sized array
    tdoas: Measured TDOA values for the selected receiver combinations
    rx: Receiver positions as an Mx2-sized array
    si: Left-hand receiver combination for tdoa computation
    sj: Right-hand receiver combination for tdoa computation

    Returns:
    np.array(N,1) with mean squared error for the given X,Y pairs
    """

    N = x.shape[0]
    M = rx.shape[0]

    # Prepare the x, y points
    P = np.vstack((x, y)).T
    P = np.reshape(P, (-1, 2, 1))
    P = np.repeat(P, M, axis=2)

    #  Prepare the Receiver matrix
    Rx = rx.T.reshape((1, 2, -1))
    Rx = np.repeat(Rx, N, axis=0)

    # Calculate distances between all points and Receivers
    drx = np.sqrt(np.sum(np.square(P - Rx), axis=1)).squeeze()
    doa = drx[:, si] - drx[:, sj]

    return np.sum(np.square(tdoas.reshape(1, -1) - doa), axis=1)


def generate_heatmap(
    tdoas: ndarray_f64,
    rx: ndarray_f64,
    xrange: Tuple[float, float],
    yrange: Tuple[float, float],
    combinations: ndarray_i64,
    step: Union[int, Tuple[int, int]] = 100,
    filter: bool = True,
    threshold: float = 0.1,
    normalize: bool = True,
) -> ndarray_f64:
    """
    Function that calculates the heatmap of the a given cost function

    Parameters:
    tdoas: np.array of shape (N,1) with TDOA values between receivers
    rx: np.array of shape (N,2) with the planar coordinates of receiver the receivers
    xrange: tuple with (xmin, xmax) values to evaluate the cost function at
    yrange: tuple with (ymin, ymax) values to evaluate the cost function at
    combinations: np.array (N,2) with sensor combinations for each of the tdoa values
    step: single float or tuple of floats with the x and y resolutions

    Returns:
    np.array(M,M) with the cost function evaluated on the grid defined by xrange, yrange &
    resolution
    """

    # Preparing the mesh to evaluate the function at
    if step is tuple:
        xnum = step[0]
        ynum = step[1]
    else:
        xnum = step
        ynum = step

    xmin, xmax = xrange
    ymin, ymax = yrange

    xvalues = np.linspace(xmin, xmax, xnum, endpoint=True)
    yvalues = np.linspace(ymin, ymax, ynum, endpoint=True)

    xmesh, ymesh = np.meshgrid(xvalues, yvalues)
    x = xmesh.flatten()
    y = ymesh.flatten()

    # We need to get the sensors involved in each combination
    si: ndarray_i64 = combinations[:, 0]
    sj: ndarray_i64 = combinations[:, 1]

    def msefun(x,y):
        mse(x, y, tdoas, rx, si, sj)

    # Core
    Z = 1 / msefun(x, y)
    if normalize:
        Z = Z / np.max(Z)

    # Remove values below threshold
    if filter:
        idx = Z > threshold
        x = x[idx]
        y = y[idx]
        Z = Z[idx]

    return np.hstack((x.reshape(-1, 1), y.reshape(-1, 1), Z.reshape(-1, 1)))

