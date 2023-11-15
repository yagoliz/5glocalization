#   Copyright (C) IMDEA Networks Institute 2023
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

from pytdoa.geodesy.geodesy import SPEED_OF_LIGHT

def dop(positions: np.ndarray, receivers: np.ndarray) -> np.ndarray:
    """
    Calculate the dilution of precision for a set of positions and receivers

    Parameters:
    positions: MxN array with M points and N coordinates
    receivers: RxN array with R receivers and N coordinates

    Returns:
    np.array(M,N) with the DOP in the N dimensions
    """

    m = positions.shape[0]
    n = positions.shape[1]
    r = receivers.shape[0]
    assert n == receivers.shape[1], "Incompatible number of dimensions for positions and receivers"
    
    # We generate a 1xnxm hypermatrix and extend it with receivers
    d_vec = np.transpose(positions.reshape(m,n,1),axes=[2,1,0]) - receivers.reshape(r,n,1)
    d = np.linalg.norm(d_vec, axis=1)
    d_vec /= (d.reshape(r,1,m)+1e-3)

    # We substract 1st row to all others
    triu_indices = np.triu_indices(r, k=1)
    H = d_vec[triu_indices[0],:,:] - d_vec[triu_indices[1],:,:]
    Q = np.linalg.inv(np.matmul(np.transpose(H,axes=[2,1,0]),np.transpose(H,axes=[2,0,1])))
    Q = np.transpose(Q,axes=[1,2,0])
    return np.hstack((Q[0,0,:].reshape(-1,1), Q[1,1,:].reshape(-1,1)))


def crlb(positions: np.ndarray, receivers: np.ndarray, c: float = SPEED_OF_LIGHT, snr: float = 10, bw: float = 100e6) -> np.ndarray:
    """
    Calculate the constant variance Cramer-Rao Lower Bound for a given set of positions and receivers

    Parameters:
    positions: MxN array with M points and N coordinates
    receivers: RxN array with R receivers and N coordinates
    c: Wave speed (by default is the speed of light)
    snr: Signal-to-Noise Ratio (in dB)
    bw: Signal bandwidth (in Hz)

    Returns:
    np.array(M,N) with the CRLB in the N dimensions for all M points
    """
    
    m = positions.shape[0]
    n = positions.shape[1]
    r = receivers.shape[0]
    assert n == receivers.shape[1], "Incompatible number of dimensions for positions and receivers"

    # We generate the variance for all sensors
    snr_w = 10 ** (snr/20)
    sigma = c**2 / (snr_w * bw**2)

    # Covariance matrix
    R = sigma*np.ones((r-1,r-1))
    di = np.diag_indices(r-1)
    R[di] += sigma
    Ri = np.linalg.inv(R)

    # We generate a 1xnxm hypermatrix and extend it with receivers
    d_vec = np.transpose(positions.reshape(m,n,1),axes=[2,1,0]) - receivers.reshape(r,n,1)
    d = np.linalg.norm(d_vec, axis=1)
    d_vec /= (d.reshape(r,1,m)+1e-3)

    mu = (d_vec[0,:,:] - d_vec[1:,:,:])
    Ri3 = np.repeat(Ri.reshape((1,r-1,r-1)),m,axis=0)
    
    J = np.matmul(np.matmul(mu.transpose([2,1,0]),Ri3),mu.transpose([2,0,1]))

    C = np.linalg.inv(J)
    Q = np.transpose(C,axes=[1,2,0])
    return np.hstack((Q[0,0,:].reshape(-1,1), Q[1,1,:].reshape(-1,1)))
