#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, matplotlib.ticker as tck, numpy as np

# create the parser
parser = ap.ArgumentParser(
    prog="Acorn Lines Plotter", description="Lines plotting script for the Quantum Acorn package.",
    formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
    add_help=False, allow_abbrev=False
)

# add the arguments for the parser
parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="This help message.")
parser.add_argument("-o", "--output", type=str, help="The output file to save the plot.")
parser.add_argument("-a", "--animate", type=int, help="Animate the plot.")

# add the plotting arguments
parser.add_argument("--colormap", type=str, nargs=2, default=["tab10", "10"], help="The colormap to use for the plot.")
parser.add_argument("--dpi", type=int, default=96, help="The DPI of the plot.")
parser.add_argument("--figsize", type=int, nargs=2, default=[6, 8], help="The dimensions of the plot in inches.")
parser.add_argument("--fontsize_label", type=int, nargs="+", help="The font size for the x and y axis labels.")
parser.add_argument("--fontsize_title", type=int, nargs="+", help="The font size for the title of each subplot.")
parser.add_argument("--fps", type=int, default=30, help="Frames per second for the animation.")
parser.add_argument("--title", type=str, nargs="+", help="The title of each subplot.")
parser.add_argument("--xlabel", type=str, nargs="+", help="The an x-axis labels for each subplot.")
parser.add_argument("--xlim", type=float, nargs="+", help="The x-axis limits for each subplot.")
parser.add_argument("--xticklabelon", type=int, nargs="+", help="Whether to show the x-axis tick labels for each subplot.")
parser.add_argument("--ylabel", type=str, nargs="+", help="The an y-axis label for each subplot.")
parser.add_argument("--ylim", type=float, nargs="+", help="The y-axis limits for each subplot.")
parser.add_argument("--yticklabelon", type=int, nargs="+", help="Whether to show the y-axis tick labels for each subplot.")

# additional file arguments
parser.add_argument("--alphas", nargs="+", help="Alpha values for the plotted lines.")
parser.add_argument("--colors", nargs="+", help="Individual colors for the plotted lines.")
parser.add_argument("--errors", nargs="+", help="Error bars for the plotted lines.")
parser.add_argument("--legends", nargs="+", help="Labels of the plotted lines.")
parser.add_argument("--markers", nargs="+", help="Markers for the plotted lines.")
parser.add_argument("--offsets", nargs="+", help="Vertical offsets for the plotted lines.")
parser.add_argument("--scales", nargs="+", help="Scales for the plotted lines.")
parser.add_argument("--styles", nargs="+", help="Styles for the plotted lines.")
parser.add_argument("--widths", nargs="+", help="Widths of the plotted lines.")
parser.add_argument("--subplots", nargs="+", help="Divide the lines into subplots.", default=["111"])

# add positional arguments and parse the arguments
parser.add_argument("files", nargs="+", help="The data files to plot."); args = parser.parse_args()

# function to load the data from the input string as well as the columns to plot and also some helper functions
def load(fname): tmp = fname.split(":"); data = np.loadtxt(tmp[0], ndmin=2, skiprows=1); return data, np.array(tmp[1].split(","), int) if len(tmp) == 2 and tmp[1] else np.arange(data.shape[1] - 1)
dcit, indc, trng = lambda d, c: ((x, j) for x, y in zip(d, c) for j in y), lambda l, e: [i for i, sl in enumerate(l) if e in sl], lambda s: range((a := list(map(int, s.split("-"))))[0], a[-1] + 1)

# load the data with plotted columns and get the total number of lines
data, data_cols = zip(*[load(file) for file in args.files]); maxnline = max(len(data_cols[i]) for i in range(len(data_cols)))

# load the error bars if provided
data_errors, cols_errors = zip(*[load(file) for file in args.errors[1:]]) if args.errors else ((), ())

# get the line indices for the error bars, colors, offsets, scales, and alphas
lind_errors  = (list(map(lambda x: int(x),  args.errors [0].split(","))) if args.errors [0] not in ["all", "every"] else np.arange(sum([len(col) for col in cols_errors]))        ) if args.errors  else []
lind_colors  = (list(map(lambda x: trng(x), args.colors [0].split(","))) if args.colors [0] not in ["all", "every"] else np.arange(len(list(dcit(data, data_cols))))[np.newaxis].T) if args.colors  else []
lind_offsets = (list(map(lambda x: trng(x), args.offsets[0].split(","))) if args.offsets[0] not in ["all", "every"] else np.arange(len(list(dcit(data, data_cols))))[np.newaxis].T) if args.offsets else []
lind_scales  = (list(map(lambda x: trng(x), args.scales [0].split(","))) if args.scales [0] not in ["all", "every"] else np.arange(len(list(dcit(data, data_cols))))[np.newaxis].T) if args.scales  else []
line_alphas  = (list(map(lambda x: trng(x), args.alphas [0].split(","))) if args.alphas [0] not in ["all", "every"] else np.arange(len(list(dcit(data, data_cols))))[np.newaxis].T) if args.alphas  else []
line_legends = (list(map(lambda x: trng(x), args.legends[0].split(","))) if args.legends[0] not in ["all", "every"] else np.arange(len(list(dcit(data, data_cols))))[np.newaxis].T) if args.legends else []
line_widths  = (list(map(lambda x: trng(x), args.widths [0].split(","))) if args.widths [0] not in ["all", "every"] else np.arange(len(list(dcit(data, data_cols))))[np.newaxis].T) if args.widths  else []
line_styles  = (list(map(lambda x: trng(x), args.styles [0].split(","))) if args.styles [0] not in ["all", "every"] else np.arange(len(list(dcit(data, data_cols))))[np.newaxis].T) if args.styles  else []
line_markers = (list(map(lambda x: trng(x), args.markers[0].split(","))) if args.markers[0] not in ["all", "every"] else np.arange(len(list(dcit(data, data_cols))))[np.newaxis].T) if args.markers else []

# add the scales and offsets to the initial lines in the data
for (f, s) in ((frame, args.animate if args.animate else 0) for frame in range((data[0].shape[1] - 1) // args.animate if args.animate else 1)):
    for i, (data_i, j) in enumerate(dcit(data, data_cols)): data_i[:, j + f * s + 1] *= float(args.scales [indc(lind_scales,  i)[0] + 1 if args.scales[0]  != "all" else 1]) if indc(lind_scales,  i) else 1
for (f, s) in ((frame, args.animate if args.animate else 0) for frame in range((data[0].shape[1] - 1) // args.animate if args.animate else 1)):
    for i, (data_i, j) in enumerate(dcit(data, data_cols)): data_i[:, j + f * s + 1] += float(args.offsets[indc(lind_offsets, i)[0] + 1 if args.offsets[0] != "all" else 1]) if indc(lind_offsets, i) else 0

# create the figure and the container for plots and error bars
fig = plt.figure(dpi=args.dpi, figsize=(args.figsize[1], args.figsize[0])); axes, plots, ebars = {}, [], []

# fill the axes dictionary
for subplot in args.subplots: axes[subplot] = (fig.add_subplot(int(subplot)) if len(subplot) == 3 else fig.add_subplot(*[int(num) for num in subplot.split("-")])) if subplot not in axes else axes[subplot]

# set the colormap
cmap = plt.get_cmap(args.colormap[0])(np.linspace(0, 1, int(args.colormap[1])))

# loop over data and its columns to plot
for i, (data_i, j) in enumerate(dcit(data, data_cols)):

    # extract the color from the colormap
    color = cmap[len([line for line in axes[args.subplots[len(plots) % len(args.subplots)]].get_lines() if not isinstance(line.get_color(), str)]) % len(cmap)]

    # set the color and alpha if provided
    alpha  = float(args.alphas [indc(line_alphas,  i)[0] + 1 if args.alphas [0] != "all" else 1]) if indc(line_alphas,  i) else 1
    color  =   str(args.colors [indc(lind_colors,  i)[0] + 1 if args.colors [0] != "all" else 1]) if indc(lind_colors,  i) else color
    legend =   str(args.legends[indc(line_legends, i)[0] + 1 if args.legends[0] != "all" else 1]) if indc(line_legends, i) else ""
    marker =   str(args.markers[indc(line_markers, i)[0] + 1 if args.markers[0] != "all" else 1]) if indc(line_markers, i) else ""
    style  =   str(args.styles [indc(line_styles,  i)[0] + 1 if args.styles [0] != "all" else 1]) if indc(line_styles,  i) else "-"
    width  = float(args.widths [indc(line_widths,  i)[0] + 1 if args.widths [0] != "all" else 1]) if indc(line_widths,  i) else 1.8

    # plot the data and append the plot to the list
    plots.append(axes[args.subplots[len(plots) % len(args.subplots)]].plot(data_i[:, 0], data_i[:, j + 1], alpha=alpha, color=color, label=legend, linestyle=style, linewidth=width, marker=marker))

# loop over the error bars and their columns
for (data_errors_i, j) in dcit(data_errors, cols_errors):

    # get the top and bottom lines
    top = plots[lind_errors[len(ebars)]][0].get_ydata() + data_errors_i[:, j + 1]
    bot = plots[lind_errors[len(ebars)]][0].get_ydata() - data_errors_i[:, j + 1]

    # extract the color
    color = plots[lind_errors[len(ebars)]][0].get_color()

    # fill the area between the top and bottom line and append the plot to the list
    ebars.append(axes[args.subplots[len(ebars) % len(args.subplots)]].fill_between(data_errors_i[:, 0], bot, top, color=color, alpha=0.2))

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

    # enable minor ticks
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator()); ax.yaxis.set_minor_locator(tck.AutoMinorLocator()); ax.tick_params(which="both", direction="in", top=True, right=True)

    # add the legend
    if any([line.get_label()[0] != "_" for line in ax.lines]): ax.legend(frameon=False)

    # set the title
    if args.title and len(args.title) > i: ax.set_title(args.title[i], fontsize=args.fontsize_title[i] if args.fontsize_title and len(args.fontsize_title) > i else 20)

    # set the axis labels
    if args.xlabel and len(args.xlabel) > i: ax.set_xlabel(args.xlabel[i], fontsize=args.fontsize_label[i] if args.fontsize_label and len(args.fontsize_label) > i else 14)
    if args.ylabel and len(args.ylabel) > i: ax.set_ylabel(args.ylabel[i], fontsize=args.fontsize_label[i] if args.fontsize_label and len(args.fontsize_label) > i else 14)

    # set the axis limits
    if args.xlim and len(args.xlim) > 2 * i and not np.isnan(args.xlim[2 * i]): ax.set_xlim(args.xlim[2 * i + 0:2 * i + 2])
    if args.ylim and len(args.ylim) > 2 * i and not np.isnan(args.ylim[2 * i]): ax.set_ylim(args.ylim[2 * i + 0:2 * i + 2])

    # enable or disable the x-axis and y-axis tick labels
    if args.xticklabelon and len(args.xticklabelon) > i: ax.xaxis.set_tick_params(labelbottom=args.xticklabelon[i])
    if args.yticklabelon and len(args.yticklabelon) > i: ax.yaxis.set_tick_params(labelleft  =args.yticklabelon[i])

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
