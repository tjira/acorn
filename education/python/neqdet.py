#!/usr/bin/env python

import argparse as ap, matplotlib.pyplot as plt, numpy as np, scipy as sp, scipy.linalg

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =======================================================================================================================================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="resmet.py", description="Numerically Exact Quantum Dynamics Educational Toolkit",
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    # add the arguments
    parser.add_argument("-w", "--guess", help="Initial guess for the wavefunction. (default: %(default)s)", type=str, default="[np.exp(-(r1-1)**2)]")
    parser.add_argument("-v", "--potential", help="Potential matrix. (default: %(default)s)", type=str, default="[[0.5*r1**2]]")
    parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)
    parser.add_argument("-i", "--iterations", help="Number of iterations. (default: %(default)s)", type=int, default=500)
    parser.add_argument("-l", "--limit", help="Distance limit of the wavefunction grid. (default: %(default)s)", type=float, default=8)
    parser.add_argument("-t", "--timestep", help="Time step of the simulation. (default: %(default)s)", type=float, default=0.1)
    parser.add_argument("-m", "--mass", help="Mass of the particle. (default: %(default)s)", type=float, default=1)
    parser.add_argument("-d", "--dimension", help="Dimensionality of the simulation. (default: %(default)s)", type=int, default=1)
    parser.add_argument("-g", "--points", help="Number of points in the wavefunction grid. (default: %(default)s)", type=int, default=128)
    parser.add_argument("-s", "--states", help="Number of states to find using imaginary time propagation. (default: %(default)s)", type=int, default=1)

    # optional flags
    parser.add_argument("-r", "--real", help="Perform real-time propagation on the last used wavefunction.", action=ap.BooleanOptionalAction)

    # parse the arguments
    args = parser.parse_args()

    # replace the variables in the string expressions
    for i in range(1, 10): args.guess     =     args.guess.replace(f"r{i}", f"r[:, {i - 1}]")
    for i in range(1, 10): args.potential = args.potential.replace(f"r{i}", f"r[:, {i - 1}]")

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # PREPARE VARIABLES FOR THE DYNAMICS ==============================================================================================================================================================

    # create the grid in real space and extract the space step
    r = np.stack(np.meshgrid(*[np.linspace(-args.limit, args.limit, args.points)] * args.dimension, indexing="ij"), axis=-1).reshape(-1, args.dimension)

    # calculate the space step and create containers for the optimized wavefunctions and time dependent observables
    dr = r[1, -1] - r[0, -1]; wfnopt, population = [], []

    # create the grid in momentum space
    k = np.stack(np.meshgrid(*[2 * np.pi * np.fft.fftfreq(args.points, dr)] * args.dimension, indexing="ij"), axis=-1).reshape(-1, args.dimension)

    # define the potential
    V = np.array(eval(args.potential)).transpose(2, 0, 1)

    # PERFORM THE PROPAGATIONS =========================================================================================================================================================================

    # iterate over the propagations
    for i in range(args.states):

        # print the propagation header
        print("%sPROPAGATION OF STATE %d " % ("\n " if i else "", i))

        # clear the time-dependent observable containers
        population.clear()

        # create the initial wavefunction from the provided guess
        psi = np.array(list(map(lambda x: x*np.ones(r.shape[0]), eval(args.guess))), dtype=complex).T

        # calculate the propagators for each point in the grid
        K = np.array([sp.linalg.expm(-0.5 * (1j if args.real else 1) * args.timestep * np.sum(k[i, :] ** 2) / args.mass * np.eye(psi.shape[1])) for i in range(r.shape[0])])
        R = np.array([sp.linalg.expm(-0.5 * (1j if args.real else 1) * args.timestep                                    * V[i]                ) for i in range(r.shape[0])])

        # print the propagation header
        print("%6s %12s %12s %12s%s" % ("ITER", "EKIN", "EPOT", "ETOT", "".join(" %12s" % ("RHO[%d, %d]" % (i, i)) for i in range(psi.shape[1]))))

        # propagate the wavefunction
        for j in range(args.iterations + 1):

            # propagate the wavefunction from second iteration
            if (j):
                psi = np.einsum("ijk,ik->ij", R, psi)
                for i in range(psi.shape[1]): psi[:, i] =  np.fft.fftn(psi[:, i].reshape(args.dimension * [args.points])).reshape(-1)
                psi = np.einsum("ijk,ik->ij", K, psi)
                for i in range(psi.shape[1]): psi[:, i] = np.fft.ifftn(psi[:, i].reshape(args.dimension * [args.points])).reshape(-1)
                psi = np.einsum("ijk,ik->ij", R, psi)

            # orthogonalize the wavefunction
            for wfn in wfnopt: psi -= np.trapz(wfn.conj() * psi, dx=dr) * wfn

            # normalize the wavefunction
            if not args.real or j == 0: psi /= np.sqrt(np.sum(psi.conj() * psi) * dr)

            # calculate the kinetic energy
            Ekin = 0; psit = psi.copy()
            for i in range(psi.shape[1]):
                psit[:, i] = np.fft.fftn(psi[:, i].reshape(args.dimension * [args.points])).reshape(-1)
            for i in range(psi.shape[1]):
                psit[:, i] = 0.5 * np.sum(k**2, axis=1) / args.mass * psit[:, i]
            for i in range(psi.shape[1]):
                psit[:, i] = np.fft.ifftn(psit[:, i].reshape(args.dimension * [args.points])).reshape(-1)
            Ekin += np.sum(psi.conj() * psit).real * dr

            # calculate the potential energy
            Epot = np.einsum("ij,ijk,ik->", psi.conj(), V, psi).real * dr
            
            # calculate the density matrix
            rho = psi.T @ psi.conj() * dr

            # assign the observables to the containers
            population.append(np.diag(rho).real)

            # print the iteration info
            if j % 100 == 0: print("%6d %12.6f %12.6f %12.6f%s" % (j, Ekin, Epot, Ekin + Epot, "".join(" %12.6f" % (rho[i, i].real) for i in range(rho.shape[0]))))

        # append the optimized wavefunction to the container
        if not args.real: wfnopt.append(psi.copy())

    # PLOT THE RESULTS ================================================================================================================================================================================

    # create the subplots
    fig, axs = plt.subplots(2, 2, figsize=(8, 6))

    # plot the population
    axs[0, 0].plot(np.arange(args.iterations + 1) * args.timestep, population)

    # set the labels
    axs[0, 0].set_xlabel("Time (a.u.)"); axs[0, 0].set_ylabel("Population")

    # show the plot
    plt.tight_layout(); plt.show()
