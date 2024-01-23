#!/usr/bin/env python

import matplotlib.animation as anm
import matplotlib.pyplot as plt

import argparse as ap
import numpy as np

import sys
import re

# define the data array and input
data, input = [], sys.stdin.read()

# define some simple useful functions
minmax = lambda df: [np.min(df), np.max(df)]
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
    energy = np.array(re.findall("E=(.*)", input), dtype=float)

    # create the figure, convert the data array and plot the potential
    [fig, ax], data = plt.subplots(), np.array(data); plt.plot(U.T[0], U.T[1])

    # calculate the Y axis minimum
    ymin, ymax = np.min(np.array([np.min(data[:, :, 1:]) + np.min(energy), np.min(U[:, 1])])), np.max(data[:, :, 1:]) + np.max(energy)

    # set the plot limits
    plt.axis((*minmax(data[np.logical_and(np.abs(data[:, :, 1]) > 1e-3, np.abs(data[:, :, 2]) > 1e-3)]), ymin - 0.02 * (ymax - ymin), ymax + 0.02 * (ymax - ymin)))

    # define the plots and animation update function
    if args.real:
        plots = [plt.plot(data[0].T[0], data[0].T[i] + energy[0])[0] for i in range(1, len(data[0][0]), 2)]
        update = lambda i: [plots[(j - 1) // 2].set_ydata(data[i].T[j] + energy[i]) for j in range(1, len(data[0][0]), 2)]
    elif args.imag:
        plots = [plt.plot(data[0].T[0], data[0].T[i] + energy[0])[0] for i in range(2, len(data[0][0]), 2)]
        update = lambda i: [plots[(j - 1) // 2].set_ydata(data[i].T[j] + energy[i]) for j in range(2, len(data[0][0]), 2)]
    elif args.dens:
        plots = [plt.plot(data[0].T[0], density(data[0].T[i], data[0].T[i + 1]) + energy[0])[0] for i in range(1, len(data[0][0]), 2)]
        update = lambda i: [plots[(j - 1) // 2].set_ydata(density(data[i].T[j], data[i].T[j + 1]) + energy[i]) for j in range(1, len(data[0][0]), 2)]
    else:
        plots = [plt.plot(data[0].T[0], data[0].T[i] + energy[0])[0] for i in range(1, len(data[0][0]))]
        update = lambda i: [plots[j - 1].set_ydata(data[i].T[j] + energy[i]) for j in range(1, len(data[0][0]))]

    # create the animation
    ani = anm.FuncAnimation(fig, update, frames=len(data), interval=30);

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
