#!/usr/bin/env python

import argparse as ap, matplotlib.pyplot as plt, numpy as np

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
    V = 0.5 * r ** 2; dr = r[1] - r[0]

    # PERFORM THE PROPAGATIONS =========================================================================================================================================================================

    # iterate over the propagations
    for i in range(args.states):

        # print the propagation header
        print("%sPROPAGATION OF STATE %d " % ("\n " if i else "", i))

        # create the initial wavefunction from the provided guess
        psi = eval(args.guess.replace("exp", "np.exp").replace("r1", "r")).astype(np.complex128);

        # calculate the propagators for each point in the grid
        K = np.exp(-0.5 * (1j if args.real else 1) * args.timestep * k ** 2 / args.mass)
        R = np.exp(-0.5 * (1j if args.real else 1) * args.timestep * V                 )

        # print the propagation header
        print("%6s %12s %12s %12s" % ("ITER", "EKIN", "EPOT", "ETOT"))

        # propagate the wavefunction
        for j in range(args.iterations + 1):

            # propagate the wavefunction from second iteration
            if (j): psi = R * np.fft.ifft(K * np.fft.fft(R * psi))

            # orthogonalize the wavefunction
            for wfn in wfnopt: psi -= np.trapz(wfn.conj() * psi, dx=dr) * wfn

            # normalize the wavefunction
            psi /= np.sqrt(np.sum(psi.conj() * psi) * dr)

            # integrate the expressions needed for kinetic and potential energy
            Ekin = np.trapz(psi.conj() * np.fft.ifft(k ** 2 * np.fft.fft(psi)), dx=dr).real
            Epot = np.trapz(psi.conj() * V * psi,                               dx=dr).real

            # calculate the total energy
            Ekin = 0.5 * Ekin / args.mass; Etot = Ekin + Epot

            # print the iteration info
            if (j % 100 == 0): print("%6d %12.6f %12.6f %12.6f" % (j, Ekin, Epot, Etot))

        # append the optimized wavefunction to the container
        wfnopt.append(psi.copy())

    # plt.plot(r, psi.real)
    # plt.plot(r, psi.imag)
    # plt.show()
