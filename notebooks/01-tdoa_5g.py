#
# We need to join the upper directory in order to access the local modules
import argparse
import os
import sys

module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)

#
import itertools
import json
import logging

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

#
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import seaborn as sns
import scipy

from scipy.ndimage import uniform_filter1d
from scipy.io import loadmat, savemat

#
from pytdoa.mlat import exact, lls, nlls
from pytdoa.pytdoa import brutefoptim, nonlinoptim, linoptim
from pytdoa.util import generate_heatmap, generate_hyperbola
from pytdoa.geodesy.geodesy import SPEED_OF_LIGHT

import pickle
from datetime import datetime


def main(N=10, type_optim="raw"):
    #
    fs = 92.16e6
    oversample = 16

    # gNodeB's
    gNB_1 = loadmat("tdoa_5g_loc/pos_gNB_1.mat")["pos_gNB_1"]
    gNB_2 = loadmat("tdoa_5g_loc/pos_gNB_2.mat")["pos_gNB_2"]
    gNB_3 = loadmat("tdoa_5g_loc/pos_gNB_3.mat")["pos_gNB_3"]

    # TOAs per gNodeB
    TOA_a = loadmat("tdoa_5g_loc/ToA_All.mat")["ToA_All"]

    # UE position
    UE_po = loadmat("tdoa_5g_loc/pos_UE_vector.mat")["pos_UE_vector"]

    cij_meas = np.array([21.64, 14.84, 6.79])

    selected_positions = np.arange(2, 17)

    # Variables to store
    cij = {}  # Dictionary per position
    error_fang = {}
    posit_fang = {}
    error_nlls = {}
    posit_nlls = {}
    error_brut = {}
    posit_brut = {}

    combination_range = itertools.combinations(np.arange(3), 2)
    combinations = np.fromiter(combination_range, dtype=np.dtype((int, 2)))

    for position in selected_positions:
        # Position reading
        gNB_1i = gNB_1[:, position]
        gNB_2i = gNB_2[:, position]
        gNB_3i = gNB_3[:, position]

        TOA_ai = TOA_a[:, position, :].squeeze().T
        TOA_ai = TOA_ai[~np.isnan(TOA_ai).any(axis=1)].T

        UE_poi = UE_po[:, position]

        # DF creation
        gNBs = pd.DataFrame(np.array([gNB_1i, gNB_2i, gNB_3i]), columns=[["x", "y"]])

        # Uniform convolution (ie moving average)
        TOA_av = uniform_filter1d(TOA_ai, size=N, mode="wrap")
        TDOA_av_total = np.array(
            [
                TOA_av[0, :] - TOA_av[1, :],
                TOA_av[0, :] - TOA_av[2, :],
                TOA_av[1, :] - TOA_av[2, :],
            ]
        )
        DDOA_av_total = TDOA_av_total / (fs * oversample) * SPEED_OF_LIGHT

        # Theoretical results
        DOA_th = np.sqrt(np.sum(np.square(UE_poi - gNBs), axis=1))
        DDOA_th = np.array(
            [DOA_th[0] - DOA_th[1], DOA_th[0] - DOA_th[2], DOA_th[1] - DOA_th[2]]
        ).reshape(-1, 1)

        cij_total = (DDOA_av_total - DDOA_th) * (fs * oversample) / SPEED_OF_LIGHT
        cij_posit = np.mean(cij_total, axis=1)

        # Storing data in dictionaries per position
        cij[f"p{position}"] = cij_total

        # Storing variables
        N_experiments = TOA_av.shape[1]
        cij_exp = np.zeros((3, N_experiments))
        error_fang_exp = np.zeros(N_experiments)
        posit_fang_exp = np.zeros((2, N_experiments))
        error_nlls_exp = np.zeros(N_experiments)
        posit_nlls_exp = np.zeros((2, N_experiments))
        error_brut_exp = np.zeros(N_experiments)
        posit_brut_exp = np.zeros((2, N_experiments))

        print(f"Position: {position} - STARTING...", end="")
        for experiment in range(N_experiments):
            Texp = TOA_av[:, experiment]

            if type_optim == "raw":
                TDOA_av = np.array(
                    [Texp[0] - Texp[1], Texp[0] - Texp[2], Texp[1] - Texp[2]]
                ).reshape(-1, 1)
            elif type_optim == "preoffset":
                TDOA_av = np.array(
                    [
                        Texp[0] - Texp[1] + cij_meas[0],
                        Texp[0] - Texp[2] + cij_meas[1],
                        Texp[1] - Texp[2] + cij_meas[2],
                    ]
                ).reshape(-1, 1)
            elif type_optim == "ideal":
                TDOA_av = np.array(
                    [
                        Texp[0] - Texp[1] - cij_posit[0],
                        Texp[0] - Texp[2] - cij_posit[1],
                        Texp[1] - Texp[2] - cij_posit[2],
                    ]
                ).reshape(-1, 1)
            DDOA_av = TDOA_av / (fs * oversample) * SPEED_OF_LIGHT

            # Exact solution
            p = linoptim(gNBs, DDOA_av[0:2], input_type="xyz").reshape(
                -1,
            )
            if len(p) == 0:
                posit_fang_exp[:, experiment] = np.nan
                error_fang_exp[experiment] = np.nan
            else:
                e = np.sqrt(np.sum((UE_poi - p) ** 2))
                posit_fang_exp[:, experiment] = p
                error_fang_exp[experiment] = e

            # NLLS solution
            p = nonlinoptim(
                gNBs,
                DDOA_av,
                combinations,
                num_dim=2,
                input_type="xyz",
                method="Newton-CG",
            ).reshape(
                -1,
            )
            e = np.sqrt(np.sum((UE_poi - p) ** 2))

            posit_nlls_exp[:, experiment] = p
            error_nlls_exp[experiment] = e

            # Brute force
            p = brutefoptim(
                gNBs,
                DDOA_av,
                combinations,
                xrange=10,
                yrange=10,
                step=0.5,
                maxiter=3,
                epsilon=1e-4,
                num_dim=2,
                input_type="xyz",
            ).reshape(
                -1,
            )
            e = np.sqrt(np.sum((UE_poi - p) ** 2))

            posit_brut_exp[:, experiment] = p
            error_brut_exp[experiment] = e

        print("FINISHED")
        error_fang[f"p{position}"] = error_fang_exp
        posit_fang[f"p{position}"] = posit_fang_exp
        error_nlls[f"p{position}"] = error_nlls_exp
        posit_nlls[f"p{position}"] = posit_nlls_exp
        error_brut[f"p{position}"] = error_brut_exp
        posit_brut[f"p{position}"] = posit_brut_exp

    results = {
        "cij": cij,
        "posit_fang": posit_fang,
        "error_fang": error_fang,
        "posit_nlls": posit_nlls,
        "error_nlls": error_nlls,
        "error_brut": error_brut,
        "posit_brut": posit_brut,
    }

    ct = int(datetime.timestamp(datetime.now()))

    with open(f"tdoa_5g_loc/results/results_{type_optim}_{N}_{ct}.pickle", "wb") as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(results, f)


if __name__ == "__main__":
    # Argparsing for or main loop
    parser = argparse.ArgumentParser(
        prog="tdoa_5g", description="Script to run TDOA localization on 5G"
    )

    parser.add_argument("-n", "--moving-average", type=int, default=10)
    parser.add_argument(
        "-t",
        "--optimization-type",
        type=str,
        choices=["raw", "preoffset", "ideal"],
        default="raw",
    )

    args = parser.parse_args()

    # Calling our main function
    N = args.moving_average
    t = args.optimization_type
    main(N=N, type_optim=t)
