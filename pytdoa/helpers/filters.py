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

def moving_average(x: np.ndarray, w: int, axis: int = 0) -> np.ndarray:
    """
    Compute the moving average of an array on a given axis

    Parameters:
    x: Array with N data points
    w: Window size
    axis: Axes on which to perform

    Returns:
    Filtered array with size N
    """

    return np.apply_along_axis(
        lambda m: np.convolve(m, np.ones(w), mode="same") / w, axis=axis, arr=x
    )


def reject_outliers(data: np.ndarray, m: int = 3, axis: int = 0) -> np.ndarray:
    """
    Reject outliers based on a simple standard deviation test

    Parameters:
    data: Array with N data points
    m: How many deviations from std to consider
    axis: Axes on which to perform

    Returns:
    Masked array with at most N data points
    """
        
    mask = ((data - np.median(data, axis=axis)) < m * np.std(data, axis=axis)).all(
        axis=1 - axis
    )
    return data[mask, :]