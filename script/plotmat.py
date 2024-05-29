#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add optional arguments
    parser.add_argument("-a", "--animate", action="store_true", help="Treat the matrix columns as time evolution with a given column step.")
    parser.add_argument("-c", "--columns", type=int, default=1, help="Number of columns that represent the function values.")
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [np.loadtxt(mat, skiprows=1) for mat in args.mats]

    # sort the matrices acccording to the first column
    mats = [mat[mat[:, 0].argsort()] for mat in mats]

    # plot the initial plots
    plots = [[plt.plot(mat[:, 0], col)[0] for col in mat[:, 1:1 + args.columns].T] for mat in mats]

    if args.animate:
        # define the update function
        update = lambda i: [[plots[j][k].set_ydata(mats[j][:, 1 + (i * args.columns + k) % (mats[j].shape[1] - 1)]) for k in range(len(plots[j]))] for j in range(len(plots))]

        # create the animation
        ani = anm.FuncAnimation(plt.gcf(), update, frames=(mats[0].shape[1] - 1) // min(args.columns, (mats[0].shape[1] - 1)), interval=30) # type: ignore

    # set the title
    plt.gcf().canvas.manager.set_window_title("Matrix Plotter") # type: ignore

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
