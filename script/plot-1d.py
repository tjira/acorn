#!/usr/bin/env python

import argparse as ap, matplotlib.animation as am, matplotlib.pyplot as pt, numpy as np

if __name__ == "__main__":
    # create the parser
    parser = ap.ArgumentParser(prog="Acorn 1D Plotter", description="One-dimensional plotting script for the Quantum Acorn package.", add_help=False)

    # add the optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-a", "--animate", type=int, help="Perform the animation with the specified column interval.")
    parser.add_argument("-l", "--legend", type=str, nargs="+", help="Add a legend to the plot.")
    parser.add_argument("-o", "--output", type=str, default="plot", help="The output file to save the plot.")
    parser.add_argument("-t", "--title", type=str, help="The title of the plot.")
    parser.add_argument("-x", "--xlabel", type=str, help="The an x-axis label.")
    parser.add_argument("-y", "--ylabel", type=str, help="The an y-axis label.")
    parser.add_argument("--png", action="store_true", help="Save the plot as a png image.")
    parser.add_argument("--pdf", action="store_true", help="Save the plot as a pdf document.")
    parser.add_argument("--mp4", action="store_true", help="Save the plot as an mp4 clip.")

    # add positional arguments
    parser.add_argument("files", nargs="+", help="The data files to plot.")

    # parse arguments and load data
    args = parser.parse_args(); data = [np.loadtxt(file.split(":")[0], skiprows=1) for file in args.files];

    # array of plotted columns for each file
    columns = [list(map(int, args.files[i].split(":")[1].split(",")) if ":" in args.files[i] else range(args.animate if args.animate else data[i].shape[1] - 1)) for i in range(len(data))]

    # calculate the limits of the plot
    xmin, xmax = min(map(lambda data: data[:, 0 ].min(), data)), max(map(lambda data: data[:, 0 ].max(), data))
    ymin, ymax = min(map(lambda data: data[:, 1:].min(), data)), max(map(lambda data: data[:, 1:].max(), data))

    # create the figure and axis and define a function that returns the line index based on the file and column
    fig, ax = pt.subplots(figsize=(8, 6)); lineind = lambda i, j: (np.cumsum(list(map(len, columns)))[i - 1] if i else 0) + j

    # initialize the lines with the initial data
    for i in range(len(data)): [ax.plot(data[i][:, 0], data[i][:, 1 + column]) for column in columns[i]]

    # plotting function
    update = lambda frame: [ax.lines[lineind(i, j)].set_ydata(data[i][:, 1 + frame * args.animate + columns[i][j]]) for j in range(len(columns[i])) for i in range(len(data))]

    # add the legend
    if args.legend: pt.legend(args.legend)

    # set the title
    if args.title: ax.set_title(args.title)

    # set the axis labels
    if args.xlabel: ax.set_xlabel(args.xlabel)
    if args.ylabel: ax.set_ylabel(args.ylabel)

    # set the limits of the plot
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin - 0.05 * (ymax - ymin), ymax + 0.05 * (ymax - ymin))

    # set the tight layout
    fig.tight_layout()

    # create the animation
    anim = am.FuncAnimation(fig, update, frames=min([(data[i].shape[1] - 1) // args.animate for i in range(len(data))]), interval=30) if args.animate else None

    # save the plot
    if args.png: fig.savefig(args.output + ".png", dpi=300)
    if args.pdf: fig.savefig(args.output + ".pdf", dpi=300)

    # save the animation
    if args.mp4: anim.save(args.output + ".mp4", writer="ffmpeg", fps=30) # type: ignore

    # show the plot
    if not args.png and not args.pdf and not args.mp4: fig.canvas.manager.set_window_title("1D Data Plotter"); pt.show(); # type: ignore
