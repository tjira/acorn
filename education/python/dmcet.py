#!/usr/bin/env python

import argparse as ap, matplotlib.pyplot as plt, numpy as np

np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%20.14f" % x))

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =======================================================================================================================================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="dmcet.py", description="Diffuse Monte Carlo Educational Toolkit",
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    # add the arguments
    parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)

    # parse the arguments
    args = parser.parse_args()

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    n = 2000;
    nmax = 5000;
    limits = [-20, 20]
    ts = 0.01
    iters = 1000

    # PERFORM THE SIMULATION ===========================================================================================================================================================================

    def V(x):
        return 0.5 * x**2

    def diffuse(walkers, ts):
        return walkers + np.sqrt(ts) * np.random.normal(size=walkers.shape)

    def branch(walkers, energy, ts):
        W = np.exp(-ts * (energy - np.mean(energy)))
        copies = list(map(lambda w: min(w, 3), (W + np.random.uniform(size=walkers.shape)).astype(int)))
        return np.repeat(walkers, copies)

    walkers, energies = np.random.uniform(limits[0], limits[1], size=(n)), []

    for step in range(iters):

        walkers = diffuse(walkers, ts)

        energy = V(walkers); eref = np.mean(energy)
        
        energies.append(eref)
        
        walkers = branch(walkers, energy, ts)

        if len(walkers) > nmax:
            walkers = walkers[np.random.choice(len(walkers), nmax, replace=False)]

        print(eref, len(walkers))
