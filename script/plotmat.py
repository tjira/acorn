#!/usr/bin/env python

import argparse as ap, matplotlib.pyplot as plt, numpy as np

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [np.loadtxt(mat, skiprows=1) for mat in args.mats]

    # plot the matrices in 1D
    [[plt.plot(mat[:, 0], col) for col in mat[:, 1:].T] for mat in mats]

    # set the title
    plt.gcf().canvas.manager.set_window_title("Matrix Plotter") # type: ignore

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
