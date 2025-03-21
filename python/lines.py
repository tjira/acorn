#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np

# create the parser
parser = ap.ArgumentParser(
    prog="Acorn Lines Plotter", description="Lines plotting script for the Quantum Acorn package.",
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
parser.add_argument("--legend", type=str, nargs="+", help="Add a legend to the drawn lines.")
parser.add_argument("--title", type=str, nargs="+", help="The title of each subplot.")
parser.add_argument("--xlabel", type=str, nargs="+", help="The an x-axis labels for each subplot.")
parser.add_argument("--xlim", type=float, nargs="+", help="The x-axis limits for each subplot.")
parser.add_argument("--ylabel", type=str, nargs="+", help="The an y-axis label for each subplot.")
parser.add_argument("--ylim", type=float, nargs="+", help="The y-axis limits for each subplot.")

# additional file arguments
parser.add_argument("--errors", nargs="+", help="Error bars for the plotted lines.")
parser.add_argument("--offsets", nargs="+", help="Vertical offsets for the plotted lines.")
parser.add_argument("--scales", nargs="+", help="Scales for the plotted lines.")
parser.add_argument("--subplots", nargs="+", help="Divide the lines into subplots.", default=["111"])

# add positional arguments and parse the arguments
parser.add_argument("files", nargs="+", help="The data files to plot."); args = parser.parse_args()

# function to load the data from the input string as well as the columns to plot
def load(fname): tmp = fname.split(":"); data = np.loadtxt(tmp[0], ndmin=2, skiprows=1); return data, np.array(tmp[1].split(","), int) if len(tmp) == 2 and tmp[1] else np.arange(data.shape[1] - 1)

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
fig = plt.figure(dpi=args.dpi, figsize=(args.figsize[1], args.figsize[0])); axes, plots, ebars = {}, [], []

# fill the axes dictionary
for subplot in args.subplots: axes[subplot] = (fig.add_subplot(int(subplot)) if len(subplot) == 3 else fig.add_subplot(*[int(num) for num in subplot.split("-")])) if subplot not in axes else axes[subplot]

# set the colormap
cmap = plt.get_cmap(args.colormap)(np.linspace(0.1, 0.9, nline) if (plt.get_cmap(args.colormap).N) > 20 else np.linspace(0, 1, plt.get_cmap(args.colormap).N))

# loop over data and its columns to plot
for i, (data_i, j) in enumerate(dcit(data, data_cols)):

    # extract the label
    label = args.legend[i] if args.legend and len(args.legend) > i else ""

    # extract the color
    color = cmap[len(axes[args.subplots[len(plots) % len(args.subplots)]].get_lines()) % len(cmap)]

    # plot the data and append the plot to the list
    plots.append(axes[args.subplots[len(plots) % len(args.subplots)]].plot(data_i[:, 0], data_i[:, j + 1], color=color, label=label))

# loop over the error bars and their columns
for (data_errors_i, j) in dcit(data_errors, cols_errors):

    # get the top and bottom lines
    top = plots[lind_errors[len(ebars)]][0].get_ydata() + data_errors_i[:, j + 1]
    bot = plots[lind_errors[len(ebars)]][0].get_ydata() - data_errors_i[:, j + 1]

    # extract the color
    color = cmap[len(axes[args.subplots[len(plots) % len(args.subplots)]].get_lines()) % len(cmap)]

    # fill the area between the top and bottom line and append the plot to the list
    ebars.append(axes[args.subplots[len(plots) % len(args.subplots)]].fill_between(data_errors_i[:, 0], bot, top, color=cmap, alpha=0.2))

# loop over the axes
for i, ax in enumerate(axes.values()):

    # list of tuples with the data and its columns to plot on the current axis
    axis_data_cols = [pair for j, pair in enumerate(dcit(data, data_cols)) if [args.subplots[k % len(args.subplots)] for k in range(len(list(dcit(data, data_cols))))][j] == list(axes.keys())[i]]

    # calculate the y limits
    ymin = np.min([data_i[:, j + frame * (args.animate if args.animate else 0) + 1].min() for frame in range((data[0].shape[1] - 1) // args.animate if args.animate else 1) for data_i, j in axis_data_cols])
    ymax = np.max([data_i[:, j + frame * (args.animate if args.animate else 0) + 1].max() for frame in range((data[0].shape[1] - 1) // args.animate if args.animate else 1) for data_i, j in axis_data_cols])

    # set the default axis limits
    ax.set_xlim(min([line.get_xdata().min() for line in ax.lines]), max([line.get_xdata().max() for line in ax.lines])); ax.set_ylim(ymin, ymax)

    # enlarge the y limit by 5 percent
    ax.set_ylim(ax.get_ylim()[0] - 0.05 * np.diff(ax.get_ylim()), ax.get_ylim()[1] + 0.05 * np.diff(ax.get_ylim()))

    # add the legend
    if any([line.get_label()[0] != "_" for line in ax.lines]): ax.legend(frameon=False)

    # set the title
    if args.title and len(args.title) > i: ax.set_title(args.title[i])

    # set the axis labels
    if args.xlabel and len(args.xlabel) > i: ax.set_xlabel(args.xlabel[i])
    if args.ylabel and len(args.ylabel) > i: ax.set_ylabel(args.ylabel[i])

    # set the axis limits
    if args.xlim and len(args.xlim) > 2 * i and not np.isnan(args.xlim[2 * i]): ax.set_xlim(args.xlim[2 * i + 0:2 * i + 2])
    if args.ylim and len(args.ylim) > 2 * i and not np.isnan(args.ylim[2 * i]): ax.set_ylim(args.ylim[2 * i + 0:2 * i + 2])

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
        ebars[i].remove(); ebars[i] = axes[args.subplots[i % len(args.subplots)]].fill_between(data_errors_i[:, 0], bot, top, color=cmap[lind_errors[i] % len(cmap)], alpha=0.2)

# create the animation
anm = anm.FuncAnimation(fig, update, frames=(data[0].shape[1] - 1) // args.animate, init_func=lambda: None, interval=1000 // args.fps) if args.animate else None

# set the window title
if not args.output: fig.canvas.manager.set_window_title("Acorn Lines Plotter")

# save or show the plot
(anm.save(args.output, fps=args.fps) if args.animate else plt.savefig(args.output)) if args.output else plt.show()
