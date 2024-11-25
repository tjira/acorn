#!/usr/bin/env python

import argparse as ap, matplotlib.animation as am, matplotlib.pyplot as pt, numpy as np

if __name__ == "__main__":
    # create the parser
    parser = ap.ArgumentParser(prog="Acorn 1D Plotter", description="One-dimensional plotting script for the Quantum Acorn package.", add_help=False)

    # add the optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-a", "--animate", type=int, help="Perform the animation with the specified column interval.")
    parser.add_argument("-b", "--bins", type=int, default=10, help="The number of bins for the histogram.")
    parser.add_argument("-d", "--dpi", type=int, default=60, help="The resolution of the plot.")
    parser.add_argument("-l", "--legend", type=str, nargs="+", help="Add a legend to the plot.")
    parser.add_argument("-o", "--output", type=str, default="plot", help="The output file to save the plot.")
    parser.add_argument("-t", "--title", type=str, help="The title of the plot.")
    parser.add_argument("-x", "--xlabel", type=str, help="The an x-axis label.")
    parser.add_argument("-y", "--ylabel", type=str, help="The an y-axis label.")
    parser.add_argument("--png", action="store_true", help="Save the plot as a png image.")
    parser.add_argument("--pdf", action="store_true", help="Save the plot as a pdf document.")
    parser.add_argument("--gif", action="store_true", help="Save the plot as a gif clip.")
    parser.add_argument("--mp4", action="store_true", help="Save the plot as an mp4 clip.")
    parser.add_argument("--histogram", action="store_true", help="Plot the histogram of the data.")

    # add positional arguments
    parser.add_argument("files", nargs="+", help="The data files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); data = [np.loadtxt(file.split(":")[0], ndmin=2, skiprows=1) for file in args.files];

    # array of plotted columns for each file
    columns = [list(map(int, args.files[i].split(":")[1].split(",")) if ":" in args.files[i] else range(args.animate if args.animate else data[i].shape[1] - (not args.histogram))) for i in range(len(data))]

    # calculate the limits of the plot
    xmin, xmax = min([data[i][:, 0                          ].min() for i in range(len(data))]) if not args.histogram else 0, max([data[i][:, 0                          ].max() for i in range(len(data))]) if not args.histogram else 0
    ymin, ymax = min([data[i][:, [j + 1 for j in columns[i]]].min() for i in range(len(data))]) if not args.histogram else 0, max([data[i][:, [j + 1 for j in columns[i]]].max() for i in range(len(data))]) if not args.histogram else 0

    # create the figure and axis and define a function that returns the line index based on the file and column
    fig, ax = pt.subplots(figsize=(8, 6)); lineind = lambda i, j: (np.cumsum(list(map(len, columns)))[i - 1] if i else 0) + j

    # initialize the lines with the initial data
    for i in (e for e in range(len(data)) if not args.histogram):
        [ax.plot(data[i][:, 0], data[i][:, 1 + column]) for column in columns[i]]

    # initialize the histogram with the initial data
    for i in (e for e in range(len(data)) if args.histogram):
        [(lambda hist, bins: ax.plot(bins[:-1], hist))(*np.histogram(data[i][:, column], bins=np.histogram(data[0][:, 0], bins=args.bins)[1], density=True)) for column in columns[i]] # type: ignore

    # plotting function
    if not args.histogram: update = lambda frame: [ax.lines[lineind(i, j)].set_ydata(data[i][:, 1 + frame * args.animate + columns[i][j]]) for j in range(len(columns[i])) for i in range(len(data))]
    if     args.histogram: update = lambda frame: [ax.lines[lineind(i, j)].set_ydata(data[i][:, 1 + frame * args.animate + columns[i][j]]) for j in range(len(columns[i])) for i in range(len(data))]

    # add the legend
    if args.legend: pt.legend(args.legend)

    # set the title
    if args.title: ax.set_title(args.title)

    # set the axis labels
    if args.xlabel: ax.set_xlabel(args.xlabel)
    if args.ylabel: ax.set_ylabel(args.ylabel)

    # failsafe against identical limits
    if xmin == xmax: xmin, xmax = xmin - 1, xmax + 1
    if ymin == ymax: ymin, ymax = ymin - 1, ymax + 1

    # set the limits of the plot
    ax.set_xlim(xmin, xmax) if not args.histogram else None; ax.set_ylim(ymin - 0.05 * (ymax - ymin), ymax + 0.05 * (ymax - ymin)) if not args.histogram else None

    # set the tight layout
    fig.tight_layout()

    # create the animation
    anim = am.FuncAnimation(fig, update, frames=min([(data[i].shape[1] - 1) // args.animate for i in range(len(data))]), interval=30) if args.animate else None

    # save the plot
    if args.png: fig.savefig(args.output + ".png", dpi=args.dpi)
    if args.pdf: fig.savefig(args.output + ".pdf", dpi=args.dpi)

    # save the animation
    if args.gif: anim.save(args.output + ".gif", writer="imagemagick", fps=30) # type: ignore
    if args.mp4: anim.save(args.output + ".mp4", writer="ffmpeg",      fps=30) # type: ignore

    # show the plot
    if not args.png and not args.pdf and not args.mp4: fig.canvas.manager.set_window_title("1D Data Plotter"); pt.show(); # type: ignore
