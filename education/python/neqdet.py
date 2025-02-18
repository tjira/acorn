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
    parser.add_argument("-g", "--guess", help="Initial guess for the wavefunction. (default: %(default)s)", type=str, default="exp(-(r1-1)**2)")
    parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)
    parser.add_argument("-i", "--iterations", help="Number of iterations. (default: %(default)s)", type=int, default=500)
    parser.add_argument("-l", "--limit", help="Distance limit of the wavefunction grid. (default: %(default)s)", type=float, default=8)
    parser.add_argument("-t", "--timestep", help="Time step of the simulation. (default: %(default)s)", type=float, default=0.1)
    parser.add_argument("-m", "--mass", help="Mass of the particle. (default: %(default)s)", type=int, default=1)
    parser.add_argument("-p", "--points", help="Number of points in the wavefunction grid. (default: %(default)s)", type=int, default=128)
    parser.add_argument("-s", "--states", help="Number of states to find using imaginary time propagation. (default: %(default)s)", type=int, default=1)

    # optional flags
    parser.add_argument("-r", "--real", help="Perform real-time propagation on the last used wavefunction.", action=ap.BooleanOptionalAction)

    # parse the arguments
    args = parser.parse_args()

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # PREPARE VARIABLES FOR THE DYNAMICS ==============================================================================================================================================================

    # container for optimized wavefunctions
    wfnopt = []

    # create the grid in real space and extract the space step
    r = np.linspace(-args.limit, args.limit, args.points)

    # create the grid in momentum space
    k = 2 * np.pi * np.fft.fftfreq(args.points, r[1] - r[0])

    # define the potential and extract the grid spacing
    # V = np.array([[0.5 * r ** 2]]).transpose(2, 0, 1); dr = r[1] - r[0]
    V = np.array([[0.01 * np.tanh(0.6 * r), 0.001 * np.exp(-r**2)], [0.001 * np.exp(-r**2), -0.01 * np.tanh(0.6 * r)]]).transpose(2, 0, 1); dr = r[1] - r[0]

    # PERFORM THE PROPAGATIONS =========================================================================================================================================================================

    # iterate over the propagations
    for i in range(args.states):

        # print the propagation header
        print("%sPROPAGATION OF STATE %d " % ("\n " if i else "", i))

        # create the initial wavefunction from the provided guess
        # psi = np.zeros((args.points, V.shape[1]), dtype=complex); psi[:, 0] = np.exp(-(r - 1)**2 + 0j * r)
        psi = np.zeros((args.points, V.shape[1]), dtype=complex); psi[:, 1] = np.exp(-(r + 10)**2 + 10j * r)

        # calculate the propagators for each point in the grid
        K = np.array([sp.linalg.expm(-0.5 * (1j if args.real else 1) * args.timestep * k[i] ** 2 / args.mass * np.eye(psi.shape[1])) for i in range(args.points)])
        R = np.array([sp.linalg.expm(-0.5 * (1j if args.real else 1) * args.timestep                         * V[i]                ) for i in range(args.points)])

        # print the propagation header
        print("%6s %12s %12s %12s%s" % ("ITER", "EKIN", "EPOT", "ETOT", "".join(" %12s" % ("RHO[%d, %d]" % (i, i)) for i in range(psi.shape[1]))))

        # propagate the wavefunction
        for j in range(args.iterations + 1):

            # propagate the wavefunction from second iteration
            if (j): psi = np.einsum("ijk,ik->ij", R, np.fft.ifft(np.einsum("ijk,ik->ij", K, np.fft.fft(np.einsum("ijk,ik->ij", R, psi), axis=0)), axis=0))

            # orthogonalize the wavefunction
            for wfn in wfnopt: psi -= np.trapz(wfn.conj() * psi, dx=dr) * wfn

            # normalize the wavefunction
            if not args.real or j == 0: psi /= np.sqrt(np.sum(psi.conj() * psi) * dr)

            # calculate the kinetic energy
            Ekin = np.sum(psi.conj() * np.fft.ifft(0.5 * k[:, None]**2 / args.mass * np.fft.fft(psi, axis=0), axis=0)).real * dr;

            # calculate the potential energy
            Epot = np.einsum("ij,ijk,ik->", psi.conj(), V, psi).real * dr
            
            # calculate the density matrix
            rho = psi.T @ psi.conj() * dr

            # print the iteration info
            if j % 100 == 0: print("%6d %12.6f %12.6f %12.6f%s" % (j, Ekin, Epot, Ekin + Epot, "".join(" %12.6f" % (rho[i, i].real) for i in range(rho.shape[0]))))

        # append the optimized wavefunction to the container
        if not args.real: wfnopt.append(psi.copy())

    # plt.plot(r, psi.real)
    # plt.plot(r, psi.imag)
    # plt.show()
