#!/usr/bin/env python

import argparse as ap, matplotlib.animation as am, matplotlib.pyplot as pt, numpy as np

if __name__ == "__main__":
    # create the parser
    parser = ap.ArgumentParser(prog="Acorn 1D Plotter", description="One-dimensional plotting script for the Quantum Acorn package.", add_help=False)

    # add the optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-a", "--animate", type=int, help="Perform the animation with the specified column interval.")
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
    columns = [np.array(list(map(int, args.files[i].split(":")[1].split(",")) if ":" in args.files[i] else [0])) for i in range(len(data))]

    # define the transform function for the data based on the columns
    transform = lambda data, columns: np.linalg.norm(data[:, columns], axis=1) if len(columns) > 1 else data[:, columns[0]]

    # calculate the limits of the plot
    vmin, vmax = min(map(lambda data: data[:, 2:].min(), data)), max(map(lambda data: data[:, 2:].max(), data))

    # create the figure and axis and define the empty plot array
    fig, ax = pt.subplots(1, len(data), figsize=(4 * len(data), 4)); ax = ax if isinstance(ax, np.ndarray) else np.array([ax]); plots = []

    # fill the plot array
    for i in range(len(data)):

        # define the independent coordinates
        x = data[i][:, 0].reshape(int(np.sqrt(data[i].shape[0])), int(np.sqrt(data[i].shape[0])))
        y = data[i][:, 1].reshape(int(np.sqrt(data[i].shape[0])), int(np.sqrt(data[i].shape[0])))

        # define the dependent coordinates
        z = transform(data[i], 2 + columns[i]).reshape(int(np.sqrt(data[i].shape[0])), int(np.sqrt(data[i].shape[0])))

        # plot the data
        plots.append(ax[i].pcolormesh(x, y, z, vmin=vmin, vmax=vmax, cmap="twilight"))

    # plotting function
    update = lambda frame: [plots[i].set_array(transform(data[i], 2 + frame * args.animate + columns[i])) for i in range(len(data))]

    # set the axis labels
    if args.xlabel: [ax.set_xlabel(args.xlabel) for ax in ax]
    if args.ylabel: [ax.set_ylabel(args.ylabel) for ax in ax]

    # set the tight layout
    fig.tight_layout()

    # create the animation
    anim = am.FuncAnimation(fig, update, frames=min([(data[i].shape[1] - 2) // args.animate for i in range(len(data))]), interval=30) if args.animate else None

    # save the plot
    if args.png: fig.savefig("plot.png", dpi=300)
    if args.pdf: fig.savefig("plot.pdf", dpi=300)

    # save the animation
    if args.mp4: anim.save("plot.mp4", writer="ffmpeg", fps=30) # type: ignore

    # show the plot
    if not args.png and not args.pdf and not args.mp4: fig.canvas.manager.set_window_title("2D Data Plotter"); pt.show(); # type: ignore
