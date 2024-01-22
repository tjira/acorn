#!/usr/bin/env python

import matplotlib.animation as anm
import matplotlib.pyplot as plt

import argparse as ap
import numpy as np

import sys
import re

# define the data arrays and input
data, densdata, input = [], [], sys.stdin.read()
xmin, xmax, ymin, ymax = 100, -100, 100, -100

# define the square of complex number
density = lambda re, im: re * re + im * im

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Acorn Wavefunction Plotter", description="Wavefunction plotting script for the Quantum Acorn package.")

    # add arguments
    parser.add_argument("--dens", action="store_true")
    parser.add_argument("--imag", action="store_true")
    parser.add_argument("--real", action="store_true")

    # parse arguments and load potential
    args = parser.parse_args(); U = np.loadtxt("U.mat")

    # parse the wavefunctions
    for line in [line for line in input.split("\n") if line]:
        data.append([]) if line[0] == "#" else data[-1].append([float(entry) for entry in line.split()])

    # extract the energy from the wavefunction file
    energy = np.array([float(re.sub("[E=,\\[\\]]", "", cell)) for cell in input.split("\n")[0].split()[len(data[0][0]) + 1:]])

    # create the figure and plot the potential
    fig, ax = plt.subplots(); plt.plot(U.T[0], U.T[1])

    # extract the wavefunction limits
    for i in range(len(data)):
        xming = np.min(np.array(data[0]).T[0][np.apply_along_axis(lambda row: row.sum(), 1, np.abs(np.array(data[i])[:, 1:])) > 0.01])
        xmaxg = np.max(np.array(data[0]).T[0][np.apply_along_axis(lambda row: row.sum(), 1, np.abs(np.array(data[i])[:, 1:])) > 0.01])
        xmin = xming if xming < xmin else xmin
        xmax = xmaxg if xmaxg > xmax else xmax

    # set the Y axis minimum
    ymin, ymax = min([np.min(np.array(data)[:, :, 1:]) + energy, np.min(U[:, 1])]), np.array(data)[:, :, 1:].max() + energy

    # set the plot limits
    plt.axis([xmin, xmax, ymin - 0.02 * (ymax - ymin), ymax + 0.02 * (ymax - ymin)])

    # define the plots and animation update function
    if args.real:
        plots = [plt.plot(np.array(data[0]).T[0], np.array(data[0]).T[i] + energy[(i - 1) // 2])[0] for i in range(1, len(data[0][0]), 2)]
        update = lambda i: [plots[(j - 1) // 2].set_ydata(np.array(data[i]).T[j] + energy[(j - 1) // 2]) for j in range(1, len(data[0][0]), 2)]
    elif args.imag:
        plots = [plt.plot(np.array(data[0]).T[0], np.array(data[0]).T[i] + energy[(i - 1) // 2])[0] for i in range(2, len(data[0][0]), 2)]
        update = lambda i: [plots[(j - 1) // 2].set_ydata(np.array(data[i]).T[j] + energy[(j - 1) // 2]) for j in range(2, len(data[0][0]), 2)]
    elif args.dens:
        plots = [plt.plot(np.array(data[0]).T[0], density(np.array(data[0]).T[i], np.array(data[0]).T[i + 1]) + energy[(i - 1) // 2])[0] for i in range(1, len(data[0][0]), 2)]
        update = lambda i: [plots[(j - 1) // 2].set_ydata(density(np.array(data[i]).T[j], np.array(data[i]).T[j + 1]) + energy[(j - 1) // 2]) for j in range(1, len(data[0][0]), 2)]
    else:
        plots = [plt.plot(np.array(data[0]).T[0], np.array(data[0]).T[i] + energy[(i - 1) // 2])[0] for i in range(1, len(data[0][0]))]
        update = lambda i: [plots[j - 1].set_ydata(np.array(data[i]).T[j] + energy[(j - 1) // 2]) for j in range(1, len(data[0][0]))]

    # create the animation
    ani = anm.FuncAnimation(fig, update, frames=len(data), interval=30);

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
