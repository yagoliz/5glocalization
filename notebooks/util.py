from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
