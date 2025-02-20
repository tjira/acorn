#!/usr/bin/env python

import argparse as ap, itertools as it, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, scipy as sp, scipy.linalg

# EXAMPLES
# ./shcdet.py -p 10 1 -r -10 0.5 -s 1 -m 2000 -t 1 -i 5000 -v "[[0.01*np.tanh(0.6*r1),0.001*np.exp(-r1**2)],[0.001*np.exp(-r1**2),-0.01*np.tanh(0.6*r1)]]" --adiabatic -n 1000 --lzsh

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =======================================================================================================================================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="resmet.py", description="Surface Hopping Classical Dynamics Educational Toolkit",
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    # add the arguments
    parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)
    parser.add_argument("-i", "--iterations", help="Number of iterations. (default: %(default)s)", type=int, default=5000)
    parser.add_argument("-m", "--mass", help="Mass of the particle. (default: %(default)s)", type=float, default=1)
    parser.add_argument("-n", "--trajectories", help="Number of classical trajectories in the simulation. (default: %(default)s)", type=int, default=1000)
    parser.add_argument("-s", "--state", help="Initial state of the trajectories. (default: %(default)s)", type=int, default=0)
    parser.add_argument("-t", "--timestep", help="Time step of the simulation. (default: %(default)s)", type=float, default=0.1)
    parser.add_argument("-p", "--momentum", help="Mean and standard deviation of momentum for each dimension. (default: %(default)s)", type=float, nargs= "+", default=[0, 1])
    parser.add_argument("-r", "--position", help="Mean and standard deviation of position for each dimension. (default: %(default)s)", type=float, nargs= "+", default=[1, 0.5])
    parser.add_argument("-u", "--adiabatic", help="Transform the results to the adiabatic basis. (default: %(default)s)", action=ap.BooleanOptionalAction)
    parser.add_argument("-v", "--potential", help="Potential matrix. (default: %(default)s)", type=str, default="[[0.5*r1**2]]")

    # surface hopping arguments
    parser.add_argument("--lzsh", help="Enable LZSH algorithm. (default: %(default)s)", action=ap.BooleanOptionalAction)

    # parse the arguments
    args = parser.parse_args()

    # replace the variables in the string expressions
    for i in range(1, 10): args.potential = args.potential.replace(f"r{i}", f"r[:, {i - 1}]")

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # PREPARE VARIABLES FOR THE DYNAMICS ==============================================================================================================================================================

    # create the initial conditions for the trajectories
    r = np.column_stack([np.random.normal(loc=mu, scale=sigma, size=args.trajectories) for mu, sigma in zip(args.position[0::2], args.position[1::2])])
    v = np.column_stack([np.random.normal(loc=mu, scale=sigma, size=args.trajectories) for mu, sigma in zip(args.momentum[0::2], args.momentum[1::2])])

    # calculate the initial velocity from mass, initial state, define the differentiation offset and some containers
    v /= args.mass; a, s = np.zeros_like(r), np.zeros((r.shape[0], args.iterations + 1), dtype=int) + args.state; diff = 0.0001; Vs, Us = [], []

    # SURFACE HOPPING ALGORITHMS ======================================================================================================================================================================

    # define the Landau-Zener Surface Hopping algorithm
    def lzsh(i, rn):

        # loop over the trajectories and states
        for (j, k) in it.product(range(args.trajectories), range(V.shape[1])):

            # calculate the first derivatives of the new state
            dk0 = (Vs[-1][j, k, k] - Vs[-2][j, k, k]) / args.timestep
            dk1 = (Vs[-2][j, k, k] - Vs[-3][j, k, k]) / args.timestep

            # calculate the first derivatives of the current state
            ds0 = (Vs[-1][j, s[j, i], s[j, i]] - Vs[-2][j, s[j, i], s[j, i]]) / args.timestep
            ds1 = (Vs[-2][j, s[j, i], s[j, i]] - Vs[-3][j, s[j, i], s[j, i]]) / args.timestep

            # calculate the second derivatives of the new and current states
            ddk = (dk0 - dk1) / args.timestep; dds = (ds0 - ds1) / args.timestep

            # check if the trajectory is in the place for a hop
            if (k == s[j, i - 1] or (dk0 - ds0) * (dk1 - ds1) > 0 or ddk * dds > 0 or ddk - dds == 0): continue

            # calculate the hopping probability
            p = np.exp(-0.5 * np.pi * np.sqrt(sqrtarg)) if ((sqrtarg := (V[j, k, k] - V[j, s[j, i], s[j, i]])**3 / (ddk - dds)) > 0) else 0

            # hop to another state if accepted
            if (not np.isnan(p) and rn < p): s[j, i:] = k

    # PERFORM THE DYNAMICS ============================================================================================================================================================================

    # print the propagation header
    print("%6s %12s %12s %12s" % ("ITER", "EKIN", "EPOT", "ETOT", ))

    # loop over the iterations
    for i in range(args.iterations + 1):

        # save the previous values
        rp = r.copy(); vp = v.copy(); ap = a.copy()

        # loop over the dimensions and propagate
        for j in (j for j in range(r.shape[1]) if i):

            # calculate offsets
            rplus = r.copy(); rplus[:, j] += diff
            rmins = r.copy(); rmins[:, j] -= diff

            # calculate the potential energy at the offsetted positions
            Vplus = np.array(eval(args.potential.replace("r", "rplus"))).transpose(2, 0, 1)
            Vmins = np.array(eval(args.potential.replace("r", "rmins"))).transpose(2, 0, 1)

            # transform the potential energy to the adiabatic basis
            if args.adiabatic: Vplus = np.einsum("ij,jk->ijk", np.linalg.eigvalsh(Vplus), np.eye(V.shape[1]))
            if args.adiabatic: Vmins = np.einsum("ij,jk->ijk", np.linalg.eigvalsh(Vmins), np.eye(V.shape[1]))

            # calculate the acceleration
            a[:, j] = 0.5 * (np.einsum("ijj->ij", Vmins)[np.arange(args.trajectories), s[:, i]] - np.einsum("ijj->ij", Vplus)[np.arange(args.trajectories), s[:, i]]) / (diff * args.mass)

            # update the positions and velocities
            v[:, j] += 0.5 * (a[:, j] + ap[:, j]) * args.timestep; r[:, j] += (v[:, j] + 0.5 * a[:, j] * args.timestep) * args.timestep

        # calculate the diabatic potential
        V = np.array(eval(args.potential)).transpose(2, 0, 1)

        # adiabatize the potential
        E, U = np.linalg.eigh(V); V = np.einsum("ij,jk->ijk", E, np.eye(V.shape[1])); Vs.append(V); Us.append(U)

        # surface hopping algorithm
        if (args.lzsh and i > 1): lzsh(i, np.random.rand())

        # calculate the potential and kinetic energy for each trajectory
        Epot = V[np.arange(args.trajectories), s[:, i], s[:, i]]; Ekin = 0.5 * args.mass * np.sum(v**2, axis=1)

        # get the mask for trajectories that hopped to another state
        mask = (s[:, i] != s[:, i - bool(i)]) & (Ekin > (Epot - V[np.arange(args.trajectories), s[:, i - bool(i)], s[:, i - bool(i)]]))

        # rescale the velocities of the trajectories that hopped to another state
        if mask.any(): v[mask] *= np.sqrt((Ekin - Epot + V[np.arange(args.trajectories), s[:, i - 1], s[:, i - 1]]) / Ekin)[mask, None]

        # recalculate the kinetic energy if the velocities were rescaled
        if mask.any(): Ekin = 0.5 * args.mass * np.sum(v**2, axis=1)

        # print the iteration
        if i % 1000 == 0: print("%6d %12.6f %12.6f %12.6f" % (i, np.mean(Ekin), np.mean(Epot), np.mean(Ekin + Epot)))

    # PRINT AND PLOT THE RESULTS ==========================================================================================================================================================================

    # print the final populations
    print("\nFINAL POPULATION: %s" % (np.bincount(s[:, -1]) / s.shape[0]))

    # create the subplots
    fig, axs = plt.subplots(1, 1, figsize=(8, 6));

    # plot the population
    axs.plot(np.arange(args.iterations + 1) * args.timestep, [np.bincount(s[:, i]) / s.shape[0] for i in range(s.shape[1])])

    # set the labels
    axs.set_xlabel("Time (a.u.)"); axs.set_ylabel("Population")

    # show the plot
    plt.tight_layout(); plt.show()
