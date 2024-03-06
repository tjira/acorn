#!/usr/bin/env python

import matplotlib.pyplot as plt

import argparse as ap
import numpy as np

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add help argument
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-d", "--dimension", type=int, default=1, help="Potential dimension.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [np.loadtxt(mat) for mat in args.mats]

    # sort by the first column
    if args.dimension == 1: mats = [mat[mat[:, 0].argsort()] for mat in mats]

    # plot the matrices
    if args.dimension == 2:
        for M in mats: plt.axes(projection="3d").plot_trisurf(M[:, 0], M[:, 1], M[:, 2])
    else:
        for M in mats: plt.plot(M[:, 0], M[:, 1])

    plt.gcf().canvas.manager.set_window_title("Matrix Plotter") # type: ignore

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
