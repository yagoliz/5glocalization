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

import logging
import numpy as np

logger = logging.getLogger("MLAT")


def fang(positions, tdoas):
    """
    Solve the TDOA Equations using Fang's exact solution
    """

    # Result variables
    x = np.empty((0))
    y = np.empty((0))

    # First part is the rotation of all the elements
    if positions.shape[0] != 3:
        raise RuntimeError("Need 3 sensor positions")

    s1 = positions[0, :].reshape(-1, 1)
    s2 = positions[1, :].reshape(-1, 1)
    s3 = positions[2, :].reshape(-1, 1)

    # Calculate the angle between s1 and s2
    theta = np.arctan2((s2[1] - s1[1]), (s2[0] - s1[0]))
    R = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    ).squeeze()

    # Rotate the vectors
    s1_rot = np.array([0.0, 0.0])

    s2_rot = (R.T @ (s2 - s1)).squeeze()
    s3_rot = (R.T @ (s3 - s1)).squeeze()

    # We extract the values for the equations
    b = s2_rot[0]
    cx = s3_rot[0]
    cy = s3_rot[1]
    c = np.sqrt(cx**2 + cy**2)

    # We extract the values for g and h
    g = ((tdoas[1] / tdoas[0]) * b - cx) / cy
    h = (
        c**2 - tdoas[1] ** 2 + tdoas[0] * tdoas[1] * (1 - (b / tdoas[0]) ** 2)
    ) / (2 * cy)

    # With this we go for the terms of the quadratic equation
    d = -(1 + g**2 - (b / tdoas[0]) ** 2)
    e = b * (1 - (b / tdoas[0]) ** 2) - 2 * g * h
    f = tdoas[0] ** 2 / 4 * (1 - (b / tdoas[0]) ** 2) ** 2 - h**2

    # Terms for x and y (positions)
    discriminant = e**2 - 4 * d * f
    if discriminant >= 0:
        xp = (-e + np.sqrt(discriminant)) / (2 * d)
        yp = g * xp + h

        xm = (-e - np.sqrt(discriminant)) / (2 * d)
        ym = g * xm + h
    else:
        return np.array([[np.inf, np.inf]])

    # Conversion to absolute coordinates
    # For the positive result
    rp = R @ np.array([xp, yp])
    rpt = rp + s1

    rpt_real = rpt.real
    rpt_norm = np.linalg.norm(rpt_real - s1) - np.linalg.norm(rpt_real - s2)

    # We need to compare whether the signs are the same for the obtained result and the observed tdoa
    if np.sign(rpt_norm) == np.sign(tdoas[0]):
        x = np.append(x, rpt_real[0])
        y = np.append(y, rpt_real[1])

    # For the negative result
    rm = R @ np.array([xm, ym])
    rmt = rm + s1

    rmt_real = rmt.real
    rmt_norm = np.linalg.norm(rmt_real - s1) - np.linalg.norm(rmt_real - s2)

    # We need to compare whether the signs are the same for the obtained result and the observed tdoa
    if np.sign(rmt_norm) == np.sign(tdoas[0]):
        if x.shape[0] != 0:
            logger.warning("Multiple solutions exist")
        x = np.append(x, rmt_real[0])
        y = np.append(y, rmt_real[1])

    return np.hstack((x.reshape(-1,1), y.reshape(-1,1)))
