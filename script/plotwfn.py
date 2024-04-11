#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, mmap, re

# define some simple useful functions
minmax = lambda df: [np.min(df), np.max(df)]
density = lambda re, im: re * re + im * im

def read(wfnfile, static=False):
    # load and map the file
    with open(wfnfile, "r+b") as file: data, filemap = [], mmap.mmap(file.fileno(), 0, prot=mmap.PROT_READ)

    # extract the header with independent variables
    info = np.array(filemap.readline().decode().split()[1:], dtype=int)
    xyz = [filemap.readline().decode().split() for _ in range(info[0])]

    # read the wavefunction values
    for _ in iter(filemap.readline, b""):
        re = np.array(filemap.readline().decode().split(), dtype=float)
        im = np.array(filemap.readline().decode().split(), dtype=float)
        data.append((re + 1j * im).reshape(*info[2:-1]))
        if static: data = data[-1:]

    # return the info, independent variables and wavefunction values
    return info, np.array(np.meshgrid(*xyz), dtype=float), np.array(data, dtype=complex)

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Wavefunction Plotter", description="Wavefunction plotting script for the Quantum Acorn package.", add_help=False)

    # add help argument
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="Show this help message and exit.")

    # add file arguments
    parser.add_argument("wfns", nargs="+", help="The wavefunction files to plot.")
    parser.add_argument("-p", "--potential", default="U.mat", help="The potential file to plot.")

    # add wfn plotting arguments
    parser.add_argument("-i", "--imag", action="store_true", help="Plot the imaginary part of the wavefunction.")
    parser.add_argument("-a", "--abs", action="store_true", help="Plot the absolute value of the wavefunction.")
    parser.add_argument("-r", "--real", action="store_true", help="Plot the real part of the wavefunction.")

    # add mode arguments
    parser.add_argument("--align", action="store_true", help="Align the wavefunctions to the potential.")
    parser.add_argument("--scale", action="store_true", help="Scale the potential to the wfn maximum.")

    # add limit arguments
    parser.add_argument("--minpt", action="store_true", help="Include the minimum of the potential.")
    parser.add_argument("--maxpt", action="store_true", help="Include the maximum of the potential.")

    # general arguments
    parser.add_argument("--static", action="store_true", help="Plot only the final frame.")

    # parse arguments, load potential and define the data
    args = parser.parse_args(); U = np.loadtxt(args.potential)

    # read the wavefunctions
    wfndata = [read(wfnfile, args.static) for wfnfile in args.wfns]

    # extract the energies
    energy = np.array([re.findall("E=(.*)", open(wfn).read()) for wfn in args.wfns], dtype=float)[:, -wfndata[0][2].shape[0]:]

    # scale the potential
    if args.scale: U.T[1:] *= np.max([np.max(np.abs(entry[2])) for entry in wfndata]) / np.max(U.T[1:])

    # align the wavefunctions to the potential
    if args.align: energy = np.repeat([U[0, 1], U[0, 2]], len(energy[0])).reshape(len(args.wfns), len(energy[0])) - energy

    if wfndata[0][0][0] == 1:
        # create the figure
        [fig, ax] = plt.subplots()

        # plot the potential
        for col in U.T[1:]: plt.plot(U.T[0], col)

        # calculate the Y axis minimum
        ymin = np.min([np.min([np.real(entry[2]), np.imag(entry[2])]) for entry in wfndata] + np.min(energy, axis=1))
        ymax = np.max([np.max([np.real(entry[2]), np.imag(entry[2])]) for entry in wfndata] + np.max(energy, axis=1))

        # include the minimum and maximum of the potential
        if args.minpt: ymin = min(ymin, np.min(U[:, 1:]))
        if args.maxpt: ymax = max(ymax, np.max(U[:, 1:]))

        # extract the effective function domain
        domain = [minmax(wfndata[i][1][0][(np.abs(wfndata[i][2]) > 1e-3).any(axis=0)]) for i in range(len(args.wfns))]

        # set the plot limits
        plt.axis((np.min(domain), np.max(domain), ymin - 0.02 * (ymax - ymin), ymax + 0.02 * (ymax - ymin)))

        # define the initial plots
        if args.real: plots = [plt.plot(wfndata[i][1][0], np.real(wfndata[i][2][-1]) + energy[i][0])[0] for i in range(len(args.wfns))]
        elif args.imag: plots = [plt.plot(wfndata[i][1][0], np.imag(wfndata[i][2][-1]) + energy[i][0])[0] for i in range(len(args.wfns))]
        elif args.abs: plots = [plt.plot(wfndata[i][1][0], np.abs(wfndata[i][2][-1]) + energy[i][0])[0] for i in range(len(args.wfns))]
        else: plots = [[plt.plot(wfndata[i][1][0], np.real(wfndata[i][2][-1]) + energy[i][0])[0], plt.plot(wfndata[i][1][0], np.imag(wfndata[i][2][0]) + energy[i][0])[0]] for i in range(len(args.wfns))]

        # define the animation update function
        def update(j):
            if args.real: [plots[i].set_ydata(np.real(wfndata[i][2][j]) + energy[i][j]) for i in range(len(plots))] #type: ignore
            elif args.imag: [plots[i].set_ydata(np.imag(wfndata[i][2][j]) + energy[i][j]) for i in range(len(plots))] #type: ignore
            elif args.abs: [plots[i].set_ydata(np.abs(wfndata[i][2][j]) + energy[i][j]) for i in range(len(plots))] #type: ignore
            else: [plots[i][0].set_ydata(np.real(wfndata[i][2][j]) + energy[i][j]) or plots[i][1].set_ydata(np.imag(wfndata[i][2][j]) + energy[i][j]) for i in range(len(plots))] #type: ignore

    elif wfndata[0][0][0] == 2:
        # create the figure
        [fig, ax] = plt.subplots((len(args.wfns) - 1) // 5 + 1, min([5, len(args.wfns)]), figsize=(3.2 * min([5, len(args.wfns)]), ((len(args.wfns) - 1) // 5 + 1) * 3))

        # define spacing and modify axes if flat
        plt.subplots_adjust(wspace=0.15); ax = ax if len(args.wfns) > 1 else np.array([ax])

        # define the plots
        plots = [ax.pcolormesh(wfndata[i][1][0], wfndata[i][1][1], np.abs(wfndata[i][2][0])) for i, ax in enumerate(ax.flat) if i < len(args.wfns)]

        # set axes limits
        [[ax.set_xlim(-6, 6), ax.set_ylim(-6, 6)] for _, ax in enumerate(ax.flat)]

        # animation update function
        def update(j):
            for i in range(len(plots)): plots[i].set_array(np.abs(wfndata[i][2][j]))

    # create the animation
    if not args.static: ani = anm.FuncAnimation(fig, update, frames=wfndata[0][0][1], interval=30); # type: ignore

    # specify the figure name
    fig.canvas.manager.set_window_title("Wavefunction Plotter") # type: ignore

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
