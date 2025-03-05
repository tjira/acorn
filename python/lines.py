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
parser.add_argument("--linewidth", type=float, default=1.8, help="The linewidth of the plot.")
parser.add_argument("--title", type=str, help="The title of the plot.")
parser.add_argument("--xlabel", type=str, help="The an x-axis label.")
parser.add_argument("--ylabel", type=str, help="The an y-axis label.")

# add positional arguments and parse the arguments
parser.add_argument("files", nargs="+", help="The data files to plot."); args = parser.parse_args()

# set the matplotlib parameters
matplotlib.rcParams["lines.linewidth"] = args.linewidth

# function to load the data from the input string as well as the columns to plot
load = lambda s: (d := np.loadtxt((p := s.split(":", 1))[0], ndmin=2, skiprows=1), np.array(p[1].split(","), int) if len(p) == 2 and p[1] else np.arange(d.shape[1] - 1) + 1)

# load the data with plotted columns and get the total number of lines
data, cols = zip(*[load(file) for file in args.files]); nline = sum(len(cols[i]) for i in range(len(cols)))

# set the colormap
cmap = plt.get_cmap(args.colormap)(np.linspace(0.1, 0.9, nline) if (ncolor := plt.get_cmap(args.colormap).N) > 20 else np.linspace(0, 1, ncolor))

# create the figure
fig, ax = plt.subplots(dpi=args.dpi, figsize=(args.figsize[1], args.figsize[0]))

# plot the data
plots = [ax.plot(data[:, 0], data[:, j], color=cmap[i]) for i, (data, j) in enumerate((data, j) for data, cols in zip(data, cols) for j in cols)]

# add the legend
if args.legend: ax.legend(args.legend, frameon=False)

# set the title
if args.title: ax.set_title(args.title)

# set the axis labels
if args.xlabel: ax.set_xlabel(args.xlabel)
if args.ylabel: ax.set_ylabel(args.ylabel)

# set the layout
fig.tight_layout()

# save or show the plot
plt.savefig(args.output) if args.output else plt.show()
