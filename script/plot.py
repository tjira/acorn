#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, pandas as pd, scipy.interpolate as si, seaborn as sns


if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add the arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-c", "--columns", type=int, default=1, help="The number of columns to plot in the matrix.")
    parser.add_argument("-d", "--dimension", type=int, default=1, help="Dimension of the data.")
    parser.add_argument("-f", "--frames", type=int, default=1, help="The number of frames to plot.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [pd.read_csv(mat, delim_whitespace=True, header=None, skiprows=1) for mat in args.mats]; data = pd.DataFrame();

    if args.dimension == 1:
        # set the column names so that the column names of function values are unique and format the data, first index is the frame, second index is the provided matrix
        for i, mat in enumerate(mats): mat.columns = [("x" if j == 0 else sum([mats[k].shape[1] - 1 for k in range(i)]) + j) for j in range(mat.shape[1])]
        data = list(map(list, zip(*[[mat[["x"] + list(mat.columns[i:args.columns + i])] for i in range(1, mat.shape[1], args.columns)] for mat in mats])))

    elif args.dimension == 2:
        # set the column names so that the column names of function values are unique and format the data, first index is the frame, second index is the provided matrix
        for i, mat in enumerate(mats): mat.columns = [("x" if j == 0 else "y" if j == 1 else sum([mats[k].shape[1] - 2 for k in range(i)]) + j - 1) for j in range(mat.shape[1])]
        data = list(map(list, zip(*[[mat[["x", "y"] + list(mat.columns[i:args.columns + i])] for i in range(2, mat.shape[1], args.columns)] for mat in mats])))

    # sort the data by first column
    data = [[frame.sort_values("x") for frame in frames] for frames in data]

    # melt all the matrices in each frame
    if args.dimension == 1: data = [pd.concat([mat.melt("x", value_name="y", var_name="var") for mat in frame]) for frame in data]

    # initialize the figure and axis
    fig = plt.figure(figsize=(8, 6)); ax = fig.add_subplot(111, projection="3d") if args.dimension == 2 else fig.add_subplot(111)

    # plot each of the matrices
    if args.dimension == 1: sns.lineplot(ax=ax, data=data[0], x="x", y="y", hue="var", palette="colorblind")
    elif args.dimension == 2:
        for mat in data[0]:
            x, y = np.meshgrid(*[np.linspace(mat[v].min(), mat[v].max(), 64) for v in ["x", "y"]])
            for i in range(2, mat.shape[1]):
                ax.plot_surface(x, y, si.griddata((mat.x, mat.y), mat.iloc[:, i], (x, y), method="cubic")) # type: ignore

    # update function for the animation
    def update(frame):
        if args.dimension == 1:
            for i in range(args.columns): ax.lines[i].set_ydata(data[frame][pd.DataFrame(data[frame])["var"] == frame * args.columns + i + 1]["y"]) # type: ignore
    
    # create the animation
    anim = anm.FuncAnimation(fig, update, frames=args.frames, interval=30) # type: ignore

    # set the plot title, tight layout, and show the plot
    fig.canvas.manager.set_window_title("Matrix Plotter"); plt.tight_layout(); plt.show() # type: ignore
