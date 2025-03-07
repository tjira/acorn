#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np

import matplotlib

# create the parser
parser = ap.ArgumentParser(
    prog="Acorn 1D Plotter", description="One-dimensional plotting script for the Quantum Acorn package.",
    formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
    add_help=False, allow_abbrev=False
)

# add the 
parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="This help message.")
parser.add_argument("-o", "--output", type=str, help="The output file to save the plot.")
parser.add_argument("-a", "--animate", type=int, help="Animate the plot.")

# add the plotting arguments
parser.add_argument("--colormap", type=str, default="tab10", help="The colormap to use for the plot.")
parser.add_argument("--dpi", type=int, default=96, help="The DPI of the plot.")
parser.add_argument("--figsize", type=int, nargs=2, default=[6, 8], help="The dimensions of the plot in inches.")
parser.add_argument("--fps", type=int, default=30, help="Frames per second for the animation.")
parser.add_argument("--legend", type=str, nargs="+", help="Add a legend to the plot.")
parser.add_argument("--title", type=str, help="The title of the plot.")
parser.add_argument("--xlabel", type=str, help="The an x-axis label.")
parser.add_argument("--xlim", type=float, nargs=2, help="The x-axis limits of the plot.")
parser.add_argument("--ylabel", type=str, help="The an y-axis label.")
parser.add_argument("--ylim", type=float, nargs=2, help="The y-axis limits of the plot.")

# additional file arguments
parser.add_argument("--errors", nargs="+", help="Error bars for the plotted lines.")
parser.add_argument("--offsets", nargs="+", help="Vertical offsets for the plotted lines.")
parser.add_argument("--scales", nargs="+", help="Scales for the plotted lines.")

# add positional arguments and parse the arguments
parser.add_argument("files", nargs="+", help="The data files to plot."); args = parser.parse_args()

# function to load the data from the input string as well as the columns to plot
load = lambda s: (d := np.loadtxt((p := s.split(":", 1))[0], ndmin=2, skiprows=1), np.array(p[1].split(","), int) if len(p) == 2 and p[1] else np.arange(d.shape[1] - 1))

# function to return the iterator for the data and columns
dcit = lambda d, c: ((d_i, j) for i, (d_i, c_i) in enumerate(zip(d, c)) for j in c_i)

# load the data with plotted columns and get the total number of lines
data, data_cols = zip(*[load(file) for file in args.files]); nline = sum(len(data_cols[i]) for i in range(len(data_cols)))

# load the error bars if provided
data_errors, cols_errors = zip(*[load(file) for file in args.errors[1:]]) if args.errors else ((), ())

# get the line indices for the error bars
lind_errors  = (list(map(lambda x: int(x), args.errors [0].split(","))) if args.errors [0] != "all" else range(sum([len(col) for col in cols_errors]))) if args.errors  else []
lind_offsets = (list(map(lambda x: int(x), args.offsets[0].split(","))) if args.offsets[0] != "all" else range(len(args.offsets) - 1                 )) if args.offsets else []
lind_scales  = (list(map(lambda x: int(x), args.scales [0].split(","))) if args.scales [0] != "all" else range(len(args.scales ) - 1                 )) if args.scales  else []

# add the scales and offsets to the initial lines in the data
for (frame, step) in ((frame, args.animate if args.animate else 0) for frame in range((data[0].shape[1] - 1) // args.animate if args.animate else 1)):
    for i, (data_i, j) in enumerate(dcit(data, data_cols)): data_i[:, j + frame * step + 1] *= float(args.scales [lind_scales .index(i) + 1]) if i in lind_scales  else 1
for (frame, step) in ((frame, args.animate if args.animate else 0) for frame in range((data[0].shape[1] - 1) // args.animate if args.animate else 1)):
    for i, (data_i, j) in enumerate(dcit(data, data_cols)): data_i[:, j + frame * step + 1] += float(args.offsets[lind_offsets.index(i) + 1]) if i in lind_offsets else 0

# create the figure and the container for plots and error bars
fig, ax = plt.subplots(dpi=args.dpi, figsize=(args.figsize[1], args.figsize[0])); plots, ebars = [], []

# set the colormap
cmap = plt.get_cmap(args.colormap)(np.linspace(0.1, 0.9, nline) if (ncolor := plt.get_cmap(args.colormap).N) > 20 else np.linspace(0, 1, ncolor))

# loop over data and its columns to plot
for data_i, j in dcit(data, data_cols):

    # plot the data and append the plot to the list
    plots.append(ax.plot(data_i[:, 0], data_i[:, j + 1], color=cmap[len(plots) % len(cmap)]))

# loop over the error bars and their columns
for (data_errors_i, j) in dcit(data_errors, cols_errors):

    # get the top and bottom lines
    top = plots[lind_errors[len(ebars)]][0].get_ydata() + data_errors_i[:, j + 1]
    bot = plots[lind_errors[len(ebars)]][0].get_ydata() - data_errors_i[:, j + 1]

    # fill the area between the top and bottom line and append the plot to the list
    ebars.append(ax.fill_between(data_errors_i[:, 0], bot, top, color=cmap[lind_errors[len(ebars)] % len(cmap)], alpha=0.2))

# calculate the y limits
ymin = np.min([[data_i[:, j + i * int(not args.animate) + 1].min() for data_i, j in dcit(data, data_cols)] for i in range((data[0].shape[1] - 1) // args.animate if args.animate else 1)])
ymax = np.max([[data_i[:, j + i * int(not args.animate) + 1].max() for data_i, j in dcit(data, data_cols)] for i in range((data[0].shape[1] - 1) // args.animate if args.animate else 1)])

# set the default axis limits
ax.set_xlim(min([data_i[:, 0 ].min() for data_i in data]), max([data_i[:, 0 ].max() for data_i in data])); ax.set_ylim(ymin, ymax)

# enlarge the y limit by 5 percent
ax.set_ylim(ax.get_ylim()[0] - 0.05 * np.diff(ax.get_ylim()), ax.get_ylim()[1] + 0.05 * np.diff(ax.get_ylim()))

# add the legend
if args.legend: ax.legend(args.legend, frameon=False)

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
        plots[i][0].set_ydata(data_i[:, j + frame * args.animate + 1])

    # loop over the error bars and their columns
    for i, (data_errors_i, j) in enumerate(dcit(data_errors, cols_errors)):

        # update the error bars
        top = plots[lind_errors[i]][0].get_ydata() + data_errors_i[:, j + 1]
        bot = plots[lind_errors[i]][0].get_ydata() - data_errors_i[:, j + 1]

        # fill the area between the top and bottom line
        ebars[i].remove(); ebars[i] = ax.fill_between(data_errors_i[:, 0], bot, top, color=cmap[lind_errors[i] % len(cmap)], alpha=0.2)

# create the animation
anm = anm.FuncAnimation(fig, update, frames=(data[0].shape[1] - 1) // args.animate, init_func=lambda: None, interval=1000 // args.fps) if args.animate else None

# save or show the plot
(anm.save(args.output, fps=args.fps) if args.animate else plt.savefig(args.output)) if args.output else plt.show()
