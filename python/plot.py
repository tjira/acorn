#!/usr/bin/env python

import argparse as ap, itertools as it, matplotlib.animation as am, matplotlib.pyplot as pt, numpy as np

def onedim(args):
    # load the data and baseline if specified
    data     = [np.loadtxt(file.split("~")[0].split(":")[0], ndmin=2, skiprows=1)                        for file in args.files]
    baseline = [np.loadtxt(file.split("~")[1].split(":")[0], ndmin=2, skiprows=1) if "~" in file else [] for file in args.files]

    # load the baseline to data map and scale the data
    basemap = [[[int(pair.split("-")[0]), int(pair.split("-")[1])] for pair in file.split("~")[1].split(":")[1].split(",")] if "~" in file else [] for file in args.files]

    # array of plotted columns for each file
    columns = [np.array(list(map(int, args.files[i].split("~")[0].split(":")[1].split(",")) if ":" in args.files[i].split("~")[0] else range(args.animate if args.animate else data[i].shape[1] - 1))) for i in range(len(data))]

    # scale the data
    for i in range(len(data)): data[i][:, 1:] *= args.scale

    # add the baseline to the data
    for (i, j, k) in ((i, j, k) for i in range(len(data)) for j in range(len(basemap[i])) for k in range(0, data[i].shape[1] - 1, args.animate if args.animate else data[i].shape[1])):
        data[i][:, 1 + k + basemap[i][j][0]] = data[i][:, 1 + k + basemap[i][j][0]] + baseline[i][:, 1 + basemap[i][j][1]] if len(baseline[i]) else 0

    index = 0
    for i in range(len(data)):
        for j in range(len(columns[i])):
            data[i][:, 1 + columns[i][j]] += args.offset[index] if len(args.offset) > index else 0; index += 1

    # get the y range to consider for the extremum calculation
    yexrange = [[1 + j + k for j in columns[i] for k in range(0, data[i].shape[1] - 1, args.animate if args.animate else data[i].shape[1])] for i in range(len(data))]

    # calculate the limits of the plot
    xmin, xmax = min([data[i][:, 0          ].min() for i in range(len(data))]), max([data[i][:, 0          ].max() for i in range(len(data))])
    ymin, ymax = min([data[i][:, yexrange[i]].min() for i in range(len(data))]), max([data[i][:, yexrange[i]].max() for i in range(len(data))])

    # create the figure and axis and define a function that returns the line index based on the file and column
    fig, ax = pt.subplots(figsize=(args.ratio[0], args.ratio[1])); lineind = lambda i, j: (np.cumsum(list(map(len, columns)))[i - 1] if i else 0) + j

    # remove the axis if requested
    if args.blank: ax.axis("off")

    # initialize the lines with the initial data
    for i in range(len(data)): [ax.plot(data[i][:, 0], data[i][:, 1 + columns[i][j]]) for j in range(len(columns[i]))]

    # plotting function
    update = lambda frame: [ax.lines[lineind(i, j)].set_ydata(data[i][:, 1 + frame * args.animate + columns[i][j]]) for j in range(len(columns[i])) for i in range(len(data))]

    # add the legend
    if args.legend: pt.legend(args.legend)

    # failsafe against identical limits
    if xmin == xmax: xmin, xmax = xmin - 1, xmax + 1
    if ymin == ymax: ymin, ymax = ymin - 1, ymax + 1

    # set the plot limits according to the arguments
    if args.domain[0] or args.domain[1]: xmin, xmax = args.domain
    if args.value[0]  or args.value[1]:  ymin, ymax = args.value

    # set the limits of the plot
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin - 0.01 * (ymax - ymin), ymax + 0.01 * (ymax - ymin))

    # create the animation
    anm = am.FuncAnimation(fig, update, frames=min([(data[i].shape[1] - 1) // args.animate for i in range(len(data))]), interval=30) if args.animate else None

    # return
    return ax, fig, anm

def twodim(args):
    # load the data
    data = [np.loadtxt(file.split("~")[0].split(":")[0], ndmin=2, skiprows=1) for file in args.files]

    # array of plotted columns for each file
    columns = [np.array(list(map(int, args.files[i].split("~")[0].split(":")[1].split(",")) if ":" in args.files[i].split("~")[0] else range(args.animate if args.animate else data[i].shape[1] - 1))) for i in range(len(data))]

    # scale the data
    for i in range(len(data)): data[i][:, 2:] *= args.scale

    # define the transform function for the data based on the columns
    transform = lambda data, columns: np.linalg.norm(data[:, columns], axis=1) if len(columns) > 1 else data[:, columns[0]]

    # calculate the limits of the plot and set the initial shape
    vmin, vmax = min(map(lambda data: data[:, 2:].min(), data)), max(map(lambda data: data[:, 2:].max(), data)); shape = [1, len(data)]

    # exctract the most appropriate shape for the data
    for i in range(2, len(data) + 1):
        if len(data) > i and len(data) % i == 0: shape = [len(data) // i, i]

    # create the figure and axis and define the empty plot array
    fig, ax = pt.subplots(shape[0], shape[1], figsize=(4 * shape[1], 4 * shape[0])); plots = []

    # reshape the axes if necessary
    ax = ax.flatten() if isinstance(ax, np.ndarray) else np.array([ax]);

    # remove the axis if requested
    if args.blank: [ax.axis("off") for ax in ax]

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

    # create the animation
    anm = am.FuncAnimation(fig, update, frames=min([(data[i].shape[1] - 2) // args.animate for i in range(len(data))]), interval=30) if args.animate else None

    # return
    return ax, fig, anm

if __name__ == "__main__":
    # create the parser
    parser = ap.ArgumentParser(prog="Acorn 1D Plotter", description="One-dimensional plotting script for the Quantum Acorn package.", add_help=False)

    # add the optional arguments
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")
    parser.add_argument("-a", "--animate", type=int, help="Perform the animation with the specified column interval.")
    parser.add_argument("-b", "--blank", action="store_true", help="Remove the axes to create a blank frame for the plot.")
    parser.add_argument("-d", "--dpi", type=int, default=100, help="The resolution of the plot.")
    parser.add_argument("-e", "--offset", type=float, nargs="+", default=[], help="Vetrical offsets for the plotted lines.")
    parser.add_argument("-f", "--fps", type=int, default=30, help="Frames per second for the animations.")
    parser.add_argument("-l", "--legend", type=str, nargs="+", help="Add a legend to the plot.")
    parser.add_argument("-m", "--domain", type=float, nargs=2, default=[0, 0], help="Domain of the plot.")
    parser.add_argument("-n", "--ndim", type=int, default=1, help="The dimension of the plot.")
    parser.add_argument("-o", "--output", type=str, default="plot", help="The output file to save the plot.")
    parser.add_argument("-r", "--ratio", type=int, nargs=2, default=[8, 6], help="The data scaling factor.")
    parser.add_argument("-s", "--scale", type=float, default=1, help="The data scaling factor.")
    parser.add_argument("-t", "--title", type=str, help="The title of the plot.")
    parser.add_argument("-v", "--value", type=float, nargs=2, default=[0, 0], help="Plotted value range.")
    parser.add_argument("-x", "--xlabel", type=str, help="The an x-axis label.")
    parser.add_argument("-y", "--ylabel", type=str, help="The an y-axis label.")
    parser.add_argument("--png", action="store_true", help="Save the plot as a png image.")
    parser.add_argument("--pdf", action="store_true", help="Save the plot as a pdf document.")
    parser.add_argument("--gif", action="store_true", help="Save the plot as a gif clip.")
    parser.add_argument("--mp4", action="store_true", help="Save the plot as an mp4 clip.")

    # add positional arguments
    parser.add_argument("files", nargs="+", help="The data files to plot.")

    # parse arguments
    args = parser.parse_args(); ax, fig, anm = onedim(args) if args.ndim == 1 else twodim(args)

    # reshape the axes if necessary
    ax = ax.flatten() if isinstance(ax, np.ndarray) else np.array([ax]); 

    # set the title
    if args.title: [ax.set_title(args.title) for ax in ax]

    # set the axis labels
    if args.xlabel: [ax.set_xlabel(args.xlabel) for ax in ax]
    if args.ylabel: [ax.set_ylabel(args.ylabel) for ax in ax]

    # set tight layout
    fig.tight_layout();

    # remove the redundant whitespace when blank is specified
    if args.blank: fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # save the plot
    if args.png: fig.savefig(args.output + ".png", dpi=args.dpi)
    if args.pdf: fig.savefig(args.output + ".pdf", dpi=args.dpi)

    # save the animation
    if args.gif: anm.save(args.output + ".gif", fps=args.fps) # type: ignore
    if args.mp4: anm.save(args.output + ".mp4", fps=args.fps) # type: ignore

    # show the plot
    if not args.png and not args.pdf and not args.gif and not args.mp4: fig.canvas.manager.set_window_title("Data Plotter"); pt.show(); # type: ignore
