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

# add the plotting arguments
parser.add_argument("--colormap", type=str, default="tab10", help="The colormap to use for the plot.")
parser.add_argument("--dpi", type=int, default=96, help="The DPI of the plot.")
parser.add_argument("--figsize", type=int, nargs=2, default=[6, 8], help="The dimensions of the plot in inches.")
parser.add_argument("--legend", type=str, nargs="+", help="Add a legend to the plot.")
parser.add_argument("--title", type=str, help="The title of the plot.")
parser.add_argument("--xlabel", type=str, help="The an x-axis label.")
parser.add_argument("--xlim", type=float, nargs=2, help="The x-axis limits of the plot.")
parser.add_argument("--ylabel", type=str, help="The an y-axis label.")
parser.add_argument("--ylim", type=float, nargs=2, help="The y-axis limits of the plot.")

# additional file arguments
parser.add_argument("-e", "--errors", nargs="+", help="Error bars for the plotted lines.")

# add positional arguments and parse the arguments
parser.add_argument("files", nargs="+", help="The data files to plot."); args = parser.parse_args()

# function to load the data from the input string as well as the columns to plot
load = lambda s: (d := np.loadtxt((p := s.split(":", 1))[0], ndmin=2, skiprows=1), np.array(p[1].split(","), int) if len(p) == 2 and p[1] else np.arange(d.shape[1] - 1))

# load the data with plotted columns and get the total number of lines
data, dcols = zip(*[load(file) for file in args.files]); nline = sum(len(dcols[i]) for i in range(len(dcols)))

# load the error bars if provided
if args.errors: errors, ecols = zip(*[load(error) for error in args.errors[1:]])

# get the line indices for the error bars
if args.errors: elines = list(map(lambda x: int(x), args.errors[0].split(","))) if args.errors[0] != "all" else np.arange(sum([len(ecol) for ecol in ecols]))

# create the figure and the container for plots and error bars
fig, ax = plt.subplots(dpi=args.dpi, figsize=(args.figsize[1], args.figsize[0])); plots, ebars = [], []

# set the colormap
cmap = plt.get_cmap(args.colormap)(np.linspace(0.1, 0.9, nline) if (ncolor := plt.get_cmap(args.colormap).N) > 20 else np.linspace(0, 1, ncolor))

# loop over data and its columns to plot
for datai, j in ((datai, j) for i, (datai, dcolsi) in enumerate(zip(data, dcols)) for j in dcolsi):

    # plot the data and append the plot to the list
    plots.append(ax.plot(datai[:, 0], datai[:, j + 1], color=cmap[len(plots) % len(cmap)]))

# loop over the error bars and their columns
for (errorsi, j) in ((errorsi, j) for i, (errorsi, ecolsi) in enumerate(zip(errors if args.errors else [], ecols if args.errors else [])) for j in ecolsi):

    # get the top and bottom lines
    top = plots[elines[len(ebars)]][0].get_ydata() + errorsi[:, j + 1]
    bot = plots[elines[len(ebars)]][0].get_ydata() - errorsi[:, j + 1]

    # fill the area between the top and bottom line and append the plot to the list
    ebars.append(ax.fill_between(errorsi[:, 0], bot, top, color=cmap[elines[len(ebars)] % len(cmap)], alpha=0.2))

# set the default axis limits
ax.set_xlim(min([datai[:, 0].min() for datai in data]), max([datai[:, 0].max() for datai in data]))

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

# save or show the plot
plt.savefig(args.output) if args.output else plt.show()
