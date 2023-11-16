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
import numpy.typing as npt

# Typing alias for better linting
ndarray_f64 = npt.NDArray[np.float64]
ndarray_i64 = npt.NDArray[np.int64]

def generate_hyperbola(
    tdoa: ndarray_f64, rx1: ndarray_f64, rx2: ndarray_f64, t: ndarray_f64
) -> ndarray_f64:
    """
    Function that calculates the hyperbola between 2 receivers given a TDOA value.

    ## Parameters:
    tdoa: TDOA value between receiver 1 and 2
    rx1: np.array of shape (1,2) with the planar coordinates of receiver 1
    rx2: np.array of shape (1,2) with the planar coordiantes of receiver 2
    t: parametric points to evaluate the hyperbola at

    ## Returns:
    np.array(2,len(t)) with the x, y coordinates of the hyperbola for a given array t
    """

    c = np.linalg.norm(rx2 - rx1) / 2

    # If estimated tdoa is larger than distance between receivers, we might be in trouble
    if np.abs(tdoa) / 2 > c:
        print(
            f"Estimated TDOA delay ({tdoa} m) is larger than distance between receivers ({c} m)"
        )
        tdoa = np.sign(tdoa) * 0.995 * c
        print("Correction TDOA delay to 0.995 RX distance")

    # Compute the hyperbola between 2 receivers
    # Note that we can calculate the canonical hyperbola and then transform and rotate
    center = (rx2 + rx1) / 2
    theta = np.arctan2((rx2[1] - rx1[1]), (rx2[0] - rx1[0]))

    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    a = tdoa / 2
    b = np.sqrt(c**2 - a**2)

    xpoints = a * np.cosh(t)
    ypoints = b * np.sinh(t)

    X_canonical = np.vstack(
        (np.append(np.flip(xpoints), xpoints), np.append(-np.flip(ypoints), ypoints))
    )
    hyp = R @ X_canonical + center.reshape((-1, 1))

    return hyp
