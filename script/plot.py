#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, pandas as pd, scipy.interpolate as si, seaborn as sns


def one(args, mats):
    # set the column names so that the column names of function values are unique and format the data, first index is the frame, second index is the provided matrix
    for i, mat in enumerate(mats): mat.columns = [("x" if j == 0 else sum([mats[k].shape[1] - 1 for k in range(i)]) + j) for j in range(mat.shape[1])]
    data = list(map(list, zip(*[[mat[["x"] + list(mat.columns[i:args.columns + i])] for i in range(1, mat.shape[1], args.columns)] for mat in mats])))

    # sort the data by first column
    data = [[frame.sort_values("x") for frame in frames] for frames in data]

    # melt the data
    data = [pd.concat([mat.melt("x", value_name="y", var_name="var") for mat in frame]) for frame in data]

    # initialize the figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # plot the lines
    sns.lineplot(ax=ax, data=data[len(data) - 1 if args.last else 0], x="x", y="y", hue="var", palette="colorblind")

    # define the animation update function
    def update(frame):
        for i in range(args.columns):
            ax.lines[i].set_ydata(data[frame][pd.DataFrame(data[frame])["var"] == frame * args.columns + i + 1]["y"])

    # set the tight layout and return
    plt.tight_layout(); return fig, update

def two(args, mats):
    # set the column names so that the column names of function values are unique and format the data, first index is the frame, second index is the provided matrix
    for i, mat in enumerate(mats): mat.columns = [("x" if j == 0 else "y" if j == 1 else sum([mats[k].shape[1] - 2 for k in range(i)]) + j - 1) for j in range(mat.shape[1])]
    data = list(map(list, zip(*[[mat[["x", "y"] + list(mat.columns[i:args.columns + i])] for i in range(2, mat.shape[1], args.columns)] for mat in mats])))

    # create the meshgrid for the data
    x, y = np.meshgrid(*[np.linspace(data[0][0][v].min(), data[0][0][v].max(), 128) for v in ["x", "y"]])

    # calculate the number of rows and columns the resulting image will have
    rows = [i for i in range(2, len(data[0]) + 1) if len(data[0]) % i == 0 and i * i <= len(data[0])]; rows = rows[-1] if len(rows) > 0 else 1; cols = len(data[0]) // rows

    # calculate min and max of the data
    zmin = np.min([[np.min(mat.iloc[:, 2]) for mat in frame] for frame in data]); zmax = np.max([[np.max(mat.iloc[:, 2]) for mat in frame] for frame in data])

    # set the heatmap parameters
    params = {"cbar":False, "xticklabels":False, "yticklabels":False, "vmin":zmin, "vmax":zmax, "rasterized":True, "cmap":"icefire"}

    # initialize the figure and axis
    fig, ax = plt.subplots(rows, cols, figsize=(3 * cols, 3 * rows))

    # plot the heatmaps
    for i, mat in enumerate(data[len(data) - 1 if args.last else 0]):
        sns.heatmap(ax=np.array([ax]).flatten()[i], data=si.griddata((mat.x, mat.y), mat.iloc[:, 2], (x, y), method="cubic"), **params)

    # define the animation update function
    def update(frame):
        for i, mat in enumerate(data[frame]): 
            ax.collections[i].remove(); sns.heatmap(ax=np.array([ax]).flatten()[i], data=si.griddata((mat.x, mat.y), mat.iloc[:, 2], (x, y), method="cubic"), **params)

    # set the tight layout and return
    plt.tight_layout(); return fig, update

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add the arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-c", "--columns", type=int, default=1, help="The number of columns to plot in the matrix.")
    parser.add_argument("-d", "--dimension", type=int, default=1, help="Dimension of the data.")
    parser.add_argument("-f", "--frames", type=int, default=1, help="The number of frames to plot.")
    parser.add_argument("--last", action="store_true", help="Display only the last frame.")
    parser.add_argument("--mp4", action="store_true", help="Save the plot as a gif.")
    parser.add_argument("--gif", action="store_true", help="Save the plot as an mp4.")
    parser.add_argument("--png", action="store_true", help="Save the plot as a png.")
    parser.add_argument("--pdf", action="store_true", help="Save the plot as a pdf.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [pd.read_csv(mat, header=None, sep="\\s+", skiprows=1) for mat in args.mats]; data = pd.DataFrame();

    # create the figure and update function
    if args.dimension == 1: fig, update = one(args, mats)
    if args.dimension == 2: fig, update = two(args, mats)

    # create the animation
    if args.frames > 1: anim = anm.FuncAnimation(fig, update, frames=args.frames, interval=30) # type: ignore

    # input/output
    if args.gif: anim.save("output.gif", writer="imagemagick", fps=30) # type: ignore
    elif args.mp4: anim.save("output.mp4", writer="ffmpeg", fps=30) # type: ignore
    elif args.png: plt.savefig("output.png", dpi=300) # type: ignore
    elif args.pdf: plt.savefig("output.pdf") # type: ignore
    else:
        fig.canvas.manager.set_window_title("Matrix Plotter"); plt.show() # type: ignore
