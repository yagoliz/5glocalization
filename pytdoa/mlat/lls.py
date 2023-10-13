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


def lls(positions, tdoas):
    """
    Solve TDOA equations using the Linear Least Squares approach
    The solutions contains the latitude and the longitude of the estimated
    transmitter position
    ---
    """

    (A, b) = getMatrices(positions, tdoas)
    result = np.linalg.lstsq(A, b)[0]
    return np.append(result[1].real, result[2].real).reshape((1,-1))


def getMatrices(positions, tdoas):
    # Initializing our dear variables
    A = np.zeros((len(tdoas), 3))
    b = np.zeros((len(tdoas), 1))

    for i in range(len(tdoas)):
        # System matrix
        A[i, 0] = -tdoas[i]
        A[i, 1] = positions[0, 0] - positions[i + 1, 0]
        A[i, 2] = positions[0, 1] - positions[i + 1, 1]

        # Solutions
        b[i] = 0.5 * (
            tdoas[i] ** 2
            + np.linalg.norm(positions[0, :]) ** 2
            - np.linalg.norm(positions[i + 1, :]) ** 2
        )

    return (A, b)
