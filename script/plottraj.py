#!/usr/bin/env python

import argparse as ap, cycler as cl, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-p", "--potential", default="U.mat", help="The potential file to plot.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [np.loadtxt(mat, ndmin=2) for mat in args.mats]

    # create the plot
    [fig, ax] = plt.subplots()

    # load and plot the columns of the potential
    U = np.loadtxt(args.potential); [ax.plot(U[:, 0], U[:, i]) for i in range(1, U.shape[1])]

    # plot the first point of the trajectory
    plots = [ax.plot(mat[0, 1], U[np.abs(U[:, 0] - mat[0, 1]).argmin(), int(mat[0, 0])], "ro")[0] for mat in mats]

    # animation update function
    def update(j):
        for i in range(len(plots)):
            plots[i].set_xdata([mats[i][j, 1]])
            plots[i].set_ydata([U[np.abs(U[:, 0] - mats[i][j, 1]).argmin(), int(mats[i][j, 0])]]) # type: ignore

    # create the animation
    ani = anm.FuncAnimation(fig, update, frames=mats[0].shape[0], interval=30); # type: ignore

    # set the title
    fig.canvas.manager.set_window_title("Trajectory Plotter") # type: ignore

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
