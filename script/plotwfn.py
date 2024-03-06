#!/usr/bin/env python

import matplotlib.animation as anm
import matplotlib.pyplot as plt

import argparse as ap
import numpy as np

import mmap, re

# define some simple useful functions
minmax = lambda df: [np.min(df), np.max(df)]
density = lambda re, im: re * re + im * im

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Wavefunction Plotter", description="Wavefunction plotting script for the Quantum Acorn package.", add_help=False)

    # add help argument
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")

    # add file arguments
    parser.add_argument("wfns", nargs="+", help="The wavefunction files to plot.")
    parser.add_argument("-d", "--dimension", type=int, default=1, help="Dimensionality of the wavefunction.")
    parser.add_argument("-p", "--potential", default="U.mat", help="The potential file to plot.")

    # add plot arguments
    parser.add_argument("--absv", action="store_true", help="Plot the wavefunction absolute value.")
    parser.add_argument("--imag", action="store_true", help="Plot the wavefunction imaginary part.")
    parser.add_argument("--real", action="store_true", help="Plot the wavefunction real part.")
    parser.add_argument("--dens", action="store_true", help="Plot the wavefunction density.")

    # add mode arguments
    parser.add_argument("--align", action="store_true", help="Align the wavefunctions to the potential.")
    parser.add_argument("--scale", action="store_true", help="Scale the potential to the wfn maximum.")

    # add limit arguments
    parser.add_argument("--minpt", action="store_true", help="Include the minimum of the potential.")
    parser.add_argument("--maxpt", action="store_true", help="Include the maximum of the potential.")

    # general arguments
    parser.add_argument("--static", action="store_true", help="Plot only the final frame.")

    # parse arguments, load potential and define the data
    args = parser.parse_args(); data, U = [], np.loadtxt(args.potential)

    # parse the wavefunctions
    for line in [line for line in sum([open(wfn).readlines() for wfn in args.wfns], []) if line]:
        data.append([]) if line[0] == "#" else data[-1].append([float(entry) for entry in line.split()])

    # extract the energy from the wavefunction file
    energy = list(re.findall("E=(.*)", "".join([open(wfn).read() for wfn in args.wfns])))

    # extract wavefunction lengths and make them equal
    lengths = [len(open(wfn).readlines()) // (len(data[0]) + 1) for wfn in args.wfns]

    # make the wavefunctions and energies equal length
    energy = [energy[min(i * max(lengths) + j, sum(lengths[:i]) + lengths[i] - 1)] for i in range(len(args.wfns)) for j in range(max(lengths))]
    data = [data[min(i * max(lengths) + j, sum(lengths[:i]) + lengths[i] - 1)] for i in range(len(args.wfns)) for j in range(max(lengths))]

    # split the data into the different wavefunctions
    data = np.array(data).reshape(len(args.wfns), len(data) // len(args.wfns), len(data[0]), len(data[0][0]))

    # split the energy into the different wavefunctions
    energy = np.array(energy, dtype=float).reshape(len(args.wfns), len(energy) // len(args.wfns))

    # scale the potential
    if args.scale: U.T[1:] *= np.max(data[:, :, :, 1:]) / np.max(U.T[1:])

    # align the wavefunctions to the potential
    if args.align: energy = np.repeat([U[0, 1], U[0, 2]], len(energy[0])).reshape(len(args.wfns), len(energy[0])) - energy

    # define the starting frame
    sf = -1 if args.static else 0

    if args.dimension == 1:
        # create the figure
        [fig, ax] = plt.subplots()

        # plot the potential
        for col in U.T[1:]: plt.plot(U.T[0], col)

        # calculate the Y axis minimum
        ymin = np.min(np.min(data[:, :, :, 1:], axis=(1, 2, 3)) + np.min(energy, axis=1))
        ymax = np.max(np.max(data[:, :, :, 1:], axis=(1, 2, 3)) + np.max(energy, axis=1))

        # include the minimum and maximum of the potential
        if args.minpt: ymin = min(ymin, np.min(U[:, 1:]))
        if args.maxpt: ymax = max(ymax, np.max(U[:, 1:]))

        # set the plot limits
        plt.axis((*minmax(data[np.logical_or(np.abs(data[:, :, :, 1]) > 1e-3, np.abs(data[:, :, :, 2]) > 1e-3)]), ymin - 0.02 * (ymax - ymin), ymax + 0.02 * (ymax - ymin)))

        # define the plots and animation update function
        if args.absv: plots = [plt.plot(data[i][0].T[0], np.sqrt(density(data[i][0].T[1], data[i][0].T[2])))[0] for i in range(len(args.wfns))]
        elif args.dens: plots = [plt.plot(data[i][0].T[0], density(data[i][0].T[1], data[i][0].T[2]))[0] for i in range(len(args.wfns))]
        elif args.imag: plots = [plt.plot(data[i][0].T[0], data[i][0].T[2])[0] for i in range(len(args.wfns))]
        elif args.real: plots = [plt.plot(data[i][0].T[0], data[i][0].T[1])[0] for i in range(len(args.wfns))]
        else: plots = [plt.plot(data[i][sf].T[0], data[i][sf].T[j + 1] + energy[i][sf])[0] for i in range(len(args.wfns)) for j in range(2)]

        # define the animation update function
        def update(i): # type: ignore
            for k in range(len(args.wfns)):
                if args.absv: plots[k].set_ydata(np.sqrt(density(data[k][i].T[1], data[k][i].T[2])) + energy[k][i])
                elif args.dens: plots[k].set_ydata(density(data[k][i].T[1], data[k][i].T[2]) + energy[k][i])
                elif args.imag: plots[k].set_ydata(data[k][i].T[2] + energy[k][i])
                elif args.real: plots[k].set_ydata(data[k][i].T[1] + energy[k][i])
                else: [plots[2 * k + j].set_ydata(data[k][i].T[j + 1] + energy[k][i]) for j in range(2)]

    elif args.dimension == 2:
        # create the figure
        [fig, ax] = plt.subplots((len(args.wfns) - 1) // 5 + 1, min([5, len(args.wfns)]), figsize=(3.2 * min([5, len(args.wfns)]), ((len(args.wfns) - 1) // 5 + 1) * 3))

        # define spacing and modify axes if flat
        plt.subplots_adjust(wspace=0.15); ax = ax if len(args.wfns) > 1 else np.array([ax])

        # get the number of gridpoints
        ngrid = np.count_nonzero(data[0][0][:, 0] == data[0][0][0, 0])

        # define the plots
        plots = [ax.pcolormesh(
            data[i][sf][:, 0].reshape(ngrid, ngrid),
            data[i][sf][:, 1].reshape(ngrid, ngrid),
            np.abs(data[i][sf][:, 2] + 1j * data[i][sf][:, 3]).reshape(ngrid, ngrid)
        ) for i, ax in enumerate(ax.flat) if i < len(args.wfns)]

        # set axes limits
        [[ax.set_xlim(-6, 6), ax.set_ylim(-6, 6)] for _, ax in enumerate(ax.flat)]

        # animation update function
        def update(j):
            for i in range(len(plots)): plots[i].set_array(np.abs(data[i][j if j < len(data[i]) else -1][:, 2] + 1j * data[i][j if j < len(data[i]) else -1][:, 3]).reshape(ngrid, ngrid))

    # create the animation
    if not args.static: ani = anm.FuncAnimation(fig, update, frames=len(data[0]), interval=30); # type: ignore

    # specify the figure name
    fig.canvas.manager.set_window_title("Wavefunction Plotter") # type: ignore

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
