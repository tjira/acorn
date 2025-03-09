#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np

# create the parser
parser = ap.ArgumentParser(
    prog="Acorn Heatmap Plotter", description="Heatmap plotting script for the Quantum Acorn package.",
    formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
    add_help=False, allow_abbrev=False
)

# add the 
parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="This help message.")
parser.add_argument("-o", "--output", type=str, help="The output file to save the plot.")
parser.add_argument("-a", "--animate", type=int, help="Animate the plot.")

# add the plotting arguments
parser.add_argument("--colormap", type=str, default="viridis", help="The colormap to use for the plot.")
parser.add_argument("--dpi", type=int, default=96, help="The DPI of the plot.")
parser.add_argument("--figsize", type=int, nargs=2, default=[6, 6], help="The dimensions of the plot in inches.")
parser.add_argument("--fps", type=int, default=30, help="Frames per second for the animation.")
parser.add_argument("--title", type=str, help="The title of the plot.")
parser.add_argument("--xlabel", type=str, help="The an x-axis label.")
parser.add_argument("--xlim", type=float, nargs=2, help="The x-axis limits of the plot.")
parser.add_argument("--ylabel", type=str, help="The an y-axis label.")
parser.add_argument("--ylim", type=float, nargs=2, help="The y-axis limits of the plot.")

# add positional arguments and parse the arguments
parser.add_argument("files", nargs="+", help="The data files to plot."); args = parser.parse_args()

# function to load the data from the input string as well as the columns to plot
load = lambda s: (d := np.loadtxt((p := s.split(":", 1))[0], ndmin=2, skiprows=1), np.array(p[1].split(","), int) if len(p) == 2 and p[1] else np.arange(d.shape[1] - 2))

# function to return the iterator for the data and columns
dcit = lambda d, c: ((d_i, j) for i, (d_i, c_i) in enumerate(zip(d, c)) for j in c_i)

# load the data with plotted columns and get the total number of lines
data, data_cols = zip(*[load(file) for file in args.files]); nline = sum(len(data_cols[i]) for i in range(len(data_cols)))

# create the figure and the container for plots
fig, ax = plt.subplots(dpi=args.dpi, figsize=(args.figsize[1], args.figsize[0])); plots = []

# calculate the y limits
zmin = np.min([[data_i[:, j + i * (args.animate if args.animate else 0) + 2].min() for data_i, j in dcit(data, data_cols)] for i in range((data[0].shape[1] - 2) // args.animate if args.animate else 1)])
zmax = np.max([[data_i[:, j + i * (args.animate if args.animate else 0) + 2].max() for data_i, j in dcit(data, data_cols)] for i in range((data[0].shape[1] - 2) // args.animate if args.animate else 1)])

# loop over data and its columns to plot
for data_i, j in dcit(data, data_cols):

    # define the independent coordinates
    x = data_i[:, 0].reshape(int(np.sqrt(data_i.shape[0])), int(np.sqrt(data_i.shape[0])))
    y = data_i[:, 1].reshape(int(np.sqrt(data_i.shape[0])), int(np.sqrt(data_i.shape[0])))

    # define the dependent coordinates
    z = data_i[:, j + 2].reshape(int(np.sqrt(data_i.shape[0])), int(np.sqrt(data_i.shape[0])))

    # plot the data and append the plot to the list
    plots.append(ax.pcolormesh(x, y, z, cmap=args.colormap, vmin=zmin, vmax=zmax))

# set the title
if args.title: ax.set_title(args.title)

# set the axis labels
if args.xlabel: ax.set_xlabel(args.xlabel)
if args.ylabel: ax.set_ylabel(args.ylabel)

# set the axis limits
if args.xlim: ax.set_xlim(args.xlim)
if args.ylim: ax.set_ylim(args.ylim)

# set the layout
fig.tight_layout()

# update function
def update(frame):

    # loop over the data and its columns to plot
    for i, (data_i, j) in enumerate(dcit(data, data_cols)):

        # update the data
        plots[i].set_array(data_i[:, j + frame * args.animate + 2])

# create the animation
anm = anm.FuncAnimation(fig, update, frames=(data[0].shape[1] - 2) // args.animate, init_func=lambda: None, interval=1000 // args.fps) if args.animate else None

# set the window title
if not args.output: fig.canvas.manager.set_window_title("Acorn Heatmap Plotter")

# save or show the plot
(anm.save(args.output, fps=args.fps) if args.animate else plt.savefig(args.output)) if args.output else plt.show()
