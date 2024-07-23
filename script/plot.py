#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.colors as mc, matplotlib.pyplot as plt, numpy as np, pandas as pd, scipy.interpolate as si, seaborn as sns


def one(args, mats):
    # set the column names of the data so that the first column is independent variable and the rest are unique dependent variables
    for i, mat in enumerate(mats): mat.columns = [("x" if j == 0 else sum([mats[k].shape[1] - 1 for k in range(i)]) + j) for j in range(mat.shape[1])]

    # split the data into frames based on the args.columns variable
    data = list(map(list, zip(*[[mat[["x"] + list(mat.columns[i:args.columns + i])] for i in range(1, mat.shape[1], args.columns)] for mat in mats])))

    # extract, sort and melt the data
    data = [pd.concat([mat[["x"] + list(mat.columns[args.extract])].sort_values("x").melt("x", value_name="y", var_name="var") for mat in frame]) for frame in data]

    # initialize the figure and axis
    fig, ax = plt.subplots(figsize=(args.resolution[0]/plt.rcParams["figure.dpi"], args.resolution[1]/plt.rcParams["figure.dpi"]));

    # set the x and y limits
    minx = min([frame.x.min() for frame in data]); maxx = max([frame.x.max() for frame in data])
    miny = min([frame.y.min() for frame in data]); maxy = max([frame.y.max() for frame in data])
    ax.set_xlim(minx, maxx); ax.set_ylim(miny - 0.05* (maxy - miny), maxy + 0.05* (maxy - miny))

    # define the color palette
    palette = sns.color_palette([args.palette], len(data[0]["var"].unique())) if args.palette in mc.BASE_COLORS | mc.CSS4_COLORS | mc.TABLEAU_COLORS | mc.XKCD_COLORS else args.palette

    # plot the lines
    sns.lineplot(ax=ax, data=data[len(data) - 1 if args.last else 0], x="x", y="y", hue="var", palette=palette, alpha=args.alpha); ax.legend(title="Column", loc="upper right")

    # remove the legend if not specified
    if not args.legend or args.image: ax.get_legend().set_visible(False)

    if args.image:
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel(""); sns.despine(ax=ax, left=True, bottom=True, right=True); plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    else: plt.tight_layout()

    # define the animation update function
    def update(frame):
        for i in range(len(args.extract)):
            ax.lines[i].set_ydata(data[frame][pd.DataFrame(data[frame])["var"] == frame * args.columns + args.extract[i]]["y"])

    # deturn the figure and update function
    return fig, update

def two(args, mats):
    # set the column names of the data so that the first two columns are independent variables and the rest are unique dependent variables
    for i, mat in enumerate(mats): mat.columns = [("x" if j == 0 else "y" if j == 1 else sum([mats[k].shape[1] - 2 for k in range(i)]) + j - 1) for j in range(mat.shape[1])]

    # split the dat into frames based on the args.columns variable
    data = list(map(list, zip(*[[mat[["x", "y"] + list(mat.columns[i:args.columns + i])] for i in range(2, mat.shape[1], args.columns)] for mat in mats])))

    # create the meshgrid for the data
    x, y = np.meshgrid(*[np.linspace(data[0][0][v].min(), data[0][0][v].max(), 128) for v in ["x", "y"]])

    # take the norm of the data if more than one column is extracted
    if len(args.extract) > 1:
        for frame in data:
            for mat in frame:
                mat.iloc[:, 2] = np.sqrt((mat.iloc[:, 2:]**2).sum(axis=1)); mat = mat.iloc[:, [0, 1, 2]]
    else: data = [[mat.iloc[:, [0, 1, 1 + args.extract[0]]] for mat in frame] for frame in data]

    # calculate the number of rows and columns the resulting image will have
    rows = [i for i in range(2, len(data[0]) + 1) if len(data[0]) % i == 0 and i * i <= len(data[0])]; rows = rows[-1] if len(rows) > 0 else 1; cols = len(data[0]) // rows

    # calculate min and max of the data
    zmin = np.min([[np.min(mat.iloc[:, 2]) for mat in frame] for frame in data]); zmax = np.max([[np.max(mat.iloc[:, 2]) for mat in frame] for frame in data])

    # set the heatmap and surface plot parameters
    hmparams = {"cbar":False, "xticklabels":False, "yticklabels":False, "vmin":zmin, "vmax":zmax, "rasterized":True, "cmap":"icefire"}
    spparams = {"rasterized":True, "cmap":"icefire", "vmin":zmin, "vmax":zmax}

    # initialize the figure and axis
    fig, ax = plt.subplots(rows, cols, figsize=(3 * cols, 3 * rows), subplot_kw={"projection": "3d" if args.surface else None})

    # set the z limits for the surface plot
    if args.surface:
        for i, axis in enumerate(np.array([ax]).flatten()): axis.set_zlim(zmin - 0.05 * (zmax - zmin), zmax + 0.05 * (zmax - zmin))

    # plot the heatmaps
    for i, mat in enumerate(data[len(data) - 1 if args.last else 0]):
        if not args.surface: sns.heatmap(ax=np.array([ax]).flatten()[i], data=si.griddata((mat.x, mat.y), mat.iloc[:, 2], (x, y), method="cubic"), **hmparams)
        else: np.array([ax]).flatten()[i].plot_trisurf(mat.x, mat.y, mat.iloc[:, 2], **spparams)

    # define the animation update function
    def update(frame):
        [[coll.remove() for coll in axis.collections] for axis in np.array([ax]).flatten()]
        for i, mat in enumerate(data[frame]): 
            if not args.surface: sns.heatmap(ax=np.array([ax]).flatten()[i], data=si.griddata((mat.x, mat.y), mat.iloc[:, 2], (x, y), method="cubic"), **hmparams)
            else: np.array([ax]).flatten()[i].plot_trisurf(mat.x, mat.y, mat.iloc[:, 2], **spparams)

    # set the tight layout and return
    plt.tight_layout(); return fig, update

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Matrix Plotter", description="Matrix plotting script for the Quantum Acorn package.", add_help=False)

    # add the arguments
    parser.add_argument("-a", "--alpha", type=float, default=1, help="The transparency of the plot.")
    parser.add_argument("-c", "--columns", type=int, default=1, help="The number of columns to plot in the matrix.")
    parser.add_argument("-d", "--dimension", type=int, default=1, help="Dimension of the data.")
    parser.add_argument("-e", "--extract", type=int, default=0, nargs="+", help="Extract the specific columns from the provided column interval.")
    parser.add_argument("-f", "--frames", type=int, default=1, help="The number of frames to plot.")
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-p", "--palette", type=str, default="colorblind", help="The color palette to use.")
    parser.add_argument("-t", "--title", type=str, default="", help="The title of the plot.")
    parser.add_argument("-r", "--resolution", type=int, nargs=2, default=[800, 600], help="The resolution of the image.")
    parser.add_argument("--image", action="store_true", help="Display only the image without frames, ticks and labels.")
    parser.add_argument("--last", action="store_true", help="Display only the last frame.")
    parser.add_argument("--legend", action="store_true", help="Display the legend.")
    parser.add_argument("--surface", action="store_true", help="Display the data as a surface plot.")
    parser.add_argument("--mp4", action="store_true", help="Save the plot as a gif.")
    parser.add_argument("--gif", action="store_true", help="Save the plot as an mp4.")
    parser.add_argument("--png", action="store_true", help="Save the plot as a png.")
    parser.add_argument("--pdf", action="store_true", help="Save the plot as a pdf.")

    # add file arguments
    parser.add_argument("mats", nargs="+", help="The matrix files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); mats = [pd.read_csv(mat, header=None, sep="\\s+", skiprows=1) for mat in args.mats]; data = pd.DataFrame();

    # set the extract variable if not provided
    if not args.extract: args.extract = range(1, args.columns + 1)

    # create the figure and update function
    if args.dimension == 1: fig, update = one(args, mats)
    if args.dimension == 2: fig, update = two(args, mats)

    # set the title of the plot
    if args.title:
        fig.gca().set_title(args.title); plt.subplots_adjust(top=0.96) # type: ignore

    # create the animation
    if args.frames > 1: anim = anm.FuncAnimation(fig, update, frames=args.frames, interval=30) # type: ignore

    # input/output
    if args.gif: anim.save("output.gif", writer="imagemagick", fps=30) # type: ignore
    elif args.mp4: anim.save("output.mp4", writer="ffmpeg", fps=30) # type: ignore
    elif args.png: plt.savefig("output.png", dpi=300) # type: ignore
    elif args.pdf: plt.savefig("output.pdf") # type: ignore
    else:
        fig.canvas.manager.set_window_title("Matrix Plotter"); plt.show() # type: ignore
