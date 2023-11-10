from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pytdoa.geodesy.geodesy import SPEED_OF_LIGHT


def make_grid(frequency: int, plot_inf: bool = False, gnbi: int = 1, gnbj: int = 2) -> Tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(2, 3, figsize=(12, 8))

    for position in range(0, 6):
        toa_measurements = pd.read_csv(
            f"comnets_data/experiments/exp{position}/{frequency}mhz.csv"
        )
        toa_measurements.drop(
            [
                "Unnamed: 0",
            ],
            axis=1,
            inplace=True,
        )

        if plot_inf:
            toa_filt = toa_measurements
        else:
            if position == 4:
                toa_measurements = toa_measurements.drop(["gNB0", "P0"], axis=1)

            toa_filt = toa_measurements.loc[
                ~((toa_measurements == float("-inf")).any(axis=1))
            ]
            toa_filt.reset_index(drop=True, inplace=True)

        i, j = position // 3, position % 3
        ax[i, j].scatter(x=np.arange(len(toa_filt)), y=toa_filt[f"gNB{gnbi}"])
        ax[i, j].scatter(x=np.arange(len(toa_filt)), y=toa_filt[f"gNB{gnbj}"])

        ax[i, j].legend([f"gNB{gnbi}", f"gNB{gnbi}"])

        ax[i, j].set_xlim([0, len(toa_filt)])

        ax[i, j].set_title(f"Position {position} - Bandwidth: {frequency} MHz")
        if i == 1:
            ax[i, j].set_xlabel("Experiment Number")
        if j == 0:
            ax[i, j].set_ylabel("Measured ToA (samples)")

    return (fig, ax)


def moving_average(x: np.ndarray, w: int, axis: int = 0) -> np.ndarray:
    return np.apply_along_axis(
        lambda m: np.convolve(m, np.ones(w), mode="same") / w, axis=axis, arr=x
    )


def reject_outliers(data: np.ndarray, m: int = 3, axis: int = 0) -> np.ndarray:
    mask = ((data - np.median(data, axis=axis)) < m * np.std(data, axis=axis)).all(
        axis=1 - axis
    )
    return data[mask, :]

def dop(positions: np.ndarray, receivers: np.ndarray) -> np.ndarray:
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







