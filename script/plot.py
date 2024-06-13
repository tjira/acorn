#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, pandas as pd, scipy.interpolate as si, seaborn as sns


def one(args, mats):
    # set the column names so that the column names of function values are unique and format the data, first index is the frame, second index is the provided matrix
    for i, mat in enumerate(mats): mat.columns = [("x" if j == 0 else sum([mats[k].shape[1] - 1 for k in range(i)]) + j) for j in range(mat.shape[1])]
    data = list(map(list, zip(*[[mat[["x"] + list(mat.columns[i:args.columns + i])] for i in range(1, mat.shape[1], args.columns)] for mat in mats])))

    # sort the data by first column
    data = [[frame.sort_values("x") for frame in frames] for frames in data]

    data = [pd.concat([mat.melt("x", value_name="y", var_name="var") for mat in frame]) for frame in data]

    # initialize the figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    sns.lineplot(ax=ax, data=data[0], x="x", y="y", hue="var", palette="colorblind")

    def update(frame):
        frame = frame + args.initial - 1
        for i in range(args.columns): ax.lines[i].set_ydata(data[frame][pd.DataFrame(data[frame])["var"] == frame * args.columns + i + 1]["y"])

    return fig, update

def two(args, mats):
    # set the column names so that the column names of function values are unique and format the data, first index is the frame, second index is the provided matrix
    for i, mat in enumerate(mats): mat.columns = [("x" if j == 0 else "y" if j == 1 else sum([mats[k].shape[1] - 2 for k in range(i)]) + j - 1) for j in range(mat.shape[1])]
    data = list(map(list, zip(*[[mat[["x", "y"] + list(mat.columns[i:args.columns + i])] for i in range(2, mat.shape[1], args.columns)] for mat in mats])))

    x, y = np.meshgrid(*[np.linspace(data[0][0][v].min(), data[0][0][v].max(), 64) for v in ["x", "y"]])

    rows = int(np.floor(np.sqrt(len(data[0])))); cols = len(data[0]) // rows; cols = cols + 1 if cols * rows < len(data[0]) else cols

    # initialize the figure and axis
    fig, ax = plt.subplots(rows, cols, figsize=(cols * 3, rows * 9 / 4))

    for i, mat in enumerate(data[args.initial - 1]): sns.heatmap(ax=ax.flatten()[i], data=si.griddata((mat.x, mat.y), mat.iloc[:, 2], (x, y), method="cubic"))

    def update(frame): pass

    return fig, update

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add the arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-c", "--columns", type=int, default=1, help="The number of columns to plot in the matrix.")
    parser.add_argument("-d", "--dimension", type=int, default=1, help="Dimension of the data.")
    parser.add_argument("-f", "--frames", type=int, default=1, help="The number of frames to plot.")
    parser.add_argument("-i", "--initial", type=int, default=1, help="Initial frame to plot.")
    parser.add_argument("--mp4", action="store_true", help="Save the plot as an mp4.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [pd.read_csv(mat, delim_whitespace=True, header=None, skiprows=1) for mat in args.mats]; data = pd.DataFrame();

    if args.dimension == 1: fig, update = one(args, mats)
    if args.dimension == 2: fig, update = two(args, mats)

    # create the animation
    anim = anm.FuncAnimation(fig, update, frames=args.frames, interval=30) # type: ignore

    # input/output
    if args.mp4: anim.save("output.mp4", writer="ffmpeg", fps=30)
    else:
        fig.canvas.manager.set_window_title("Matrix Plotter"); plt.tight_layout(); plt.show() # type: ignore
