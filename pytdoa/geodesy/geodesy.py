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

from dataclasses import dataclass

import numpy as np

# Speed of light
SPEED_OF_LIGHT = 299792458

# Geodetic class
@dataclass
class geoC:
    # Basic Earth model constants
    WGS84_A: float = 6378137.0
    WGS84_F: float = 1.0 / 298.257223563
    WGS84_B: float = WGS84_A * (1 - WGS84_F)
    WGS84_ECC_SQ: float = 1 - WGS84_B * WGS84_B / (WGS84_A * WGS84_A)
    WGS84_ECC: float = np.sqrt(WGS84_ECC_SQ)

    # Average radius for a spherical Earth
    SPHERICAL_R: float = 6371e3
    CIRCUMFERENCE: float = 2 * np.pi * SPHERICAL_R

    # Some derived values
    wgs84_ep: float = np.sqrt((WGS84_A**2 - WGS84_B**2) / WGS84_B**2)
    wgs84_ep2_b: float = wgs84_ep**2 * WGS84_B
    wgs84_e2_a: float = WGS84_ECC_SQ * WGS84_A


def dist3fromllh(llh0, llh1):
    """
    Computing 3D distance from latitude, longitude and altitude

    Parameters:
    llh0: numpy array of shape(N,3)
    llh1: numpy array of shape(N,3)

    Returns:
    np.array(N,3): 3D distance
    """

    p0 = llh2ecef(llh0)
    p1 = llh2ecef(llh1)

    return ecef_distance(p0, p1)


def ecef_distance(p0, p1):
    """
    Distance from ECEF/Euclidean (XYZ) coordinates

    Parameters:
    p0: numpy array of shape(N,num_dim)
    p1: numpy array of shape(N,num_dim)

    Returns:
    np.array(N,num_dim): 2D,3D distance
    """

    return np.sqrt(np.sum((p0 - p1) ** 2, axis=1))


def llh2ecef(llh):
    """
    Compute ECEF distance from latitude, longitude and altitude

    Parameters:
    llh: numpy array of shape(N,3)

    Returns:
    np.array(N,3): ECEF distance
    """

    lat = deg2rad(llh[:, 0])
    lon = deg2rad(llh[:, 1])
    alt = llh[:, 2]

    slat = np.sin(lat)
    slon = np.sin(lon)
    clat = np.cos(lat)
    clon = np.cos(lon)

    d = np.sqrt(1 - (slat * slat * geoC.WGS84_ECC_SQ))
    rn = geoC.WGS84_A / d

    x = (rn + alt) * clat * clon
    y = (rn + alt) * clat * slon
    z = (rn * (1 - geoC.WGS84_ECC_SQ) + alt) * slat

    return np.vstack((x, y, z)).T


def ecef2llh(ecef):
    """
    Compute latitude, longitude and altitude from ECEF distance

    Parameters:
    ecef: numpy array of shape(N,3)

    Returns:
    np.array(N,3): latitude, longitude, altitude
    """

    x = ecef[:, 0]
    y = ecef[:, 1]
    z = ecef[:, 2]

    lon = np.arctan2(y, x)

    p = np.sqrt(x**2 + y**2)
    th = np.arctan2(geoC.WGS84_A * z, geoC.WGS84_B * p)
    lat = np.arctan2(
        z + geoC.wgs84_ep2_b * np.sin(th) ** 3, p - geoC.wgs84_e2_a * np.cos(th) ** 3
    )

    N = geoC.WGS84_A / np.sqrt(1 - geoC.WGS84_ECC_SQ * np.sin(lat) ** 2)
    alt = p / np.cos(lat) - N

    return np.vstack((rad2deg(lat), rad2deg(lon), alt)).T


def rad2deg(radian):
    return radian * 180.0 / np.pi


def deg2rad(angle):
    return angle * np.pi / 180.0


def havdist(ll0, ll1):
    lat0 = deg2rad(ll0[:, 0])
    lon0 = deg2rad(ll0[:, 1])

    lat1 = deg2rad(ll1[:, 0])
    lon1 = deg2rad(ll1[:, 1])

    hav = haversine(lat1 - lat0) + np.cos(lat0) * np.cos(lat1) * haversine(lon1 - lon0)
    return 2 * geoC.WGS84_A * np.arcsin(np.sqrt(hav))


def haversine(theta):
    return np.sin(theta / 2) ** 2


def latlon2xy(lat, lon, ref_lat, ref_lon):
    y = (lat - ref_lat) / 360 * geoC.CIRCUMFERENCE
    x = (lon - ref_lon) / 360 * np.cos(ref_lat * np.pi / 180) * geoC.CIRCUMFERENCE

    return np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))


def xy2latlon(x, y, ref_lat, ref_lon):
    lat = (y * 360 / geoC.CIRCUMFERENCE) + ref_lat
    lon = ((x * 360) / (geoC.CIRCUMFERENCE * np.cos(ref_lat * np.pi / 180))) + ref_lon

    return np.hstack((lat.reshape(-1, 1), lon.reshape(-1, 1)))
