#!/usr/bin/env python

import argparse as ap, cycler as cl, matplotlib.pyplot as plt, numpy as np

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-d", "--dimension", type=int, default=1, help="Potential dimension.")
    parser.add_argument("-n", "--normalize", action="store_true", help="Normalize provided data.")
    parser.add_argument("-c", "--complex", action="store_true", help="Treat the data as complex numbers and plot their absolute value.")
    parser.add_argument("--colors", nargs="+", help="Colors for the data.")
    parser.add_argument("--legend", nargs="+", help="Legend for the data.")
    parser.add_argument("--xlim", nargs=2, type=float, help="Limits on the horizontal axis.")
    parser.add_argument("--ylim", nargs=2, type=float, help="Limits on the vertical axis.")
    parser.add_argument("--xlab", help="Label for the horizontal axis.")
    parser.add_argument("--ylab", help="Label for the vertical axis.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [np.loadtxt(mat) for mat in args.mats]

    # normalize the data
    for i in (i for i in range(len(mats)) if args.normalize): mats[i][:, 1] /= np.max(mats[i][:, 1])

    # sort by the first column
    if args.dimension == 1: mats = [mat[mat[:, 0].argsort()] for mat in mats]

    # add colors if provided
    if args.colors: plt.gca().set_prop_cycle(cl.cycler(color=args.colors))

    # convert matrices to complex if required
    if args.complex: mats = [np.hstack((mat[:, 0].reshape(-1, 1), mat[:, 1::2] + 1j * mat[:, 2::2])) for mat in mats]

    # plot the matrices in 1D
    [[plt.plot(np.real(M[:, 0]), np.abs(C) if args.complex else C) for C in M[:, 1:].T] for M in mats if args.dimension == 1]

    # plot the matrices in 2D
    [plt.axes(projection="3d").plot_trisurf(M[:, 0], M[:, 1], M[:, 2]) for M in mats if args.dimension == 2]

    # add legend, labels and limits if provided
    if args.legend: plt.legend(args.legend)
    if args.xlab: plt.xlabel(args.xlab)
    if args.ylab: plt.ylabel(args.ylab)
    if args.xlim: plt.xlim(args.xlim)
    if args.ylim: plt.ylim(args.ylim)

    # set the title
    plt.gcf().canvas.manager.set_window_title("Matrix Plotter") # type: ignore

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
