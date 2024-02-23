#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np

def energy(wfn):
    Ek = 0.5 * np.conj(wfn) * np.fft.ifft(k**2 * np.fft.fft(wfn))
    Er = np.conj(wfn) * V * wfn; return np.sum(Ek + Er).real * dx

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =======================================================================================================================================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="eqdet.py", description="Exact Quantum Adiabatic Dynamics Educational Toolkit", add_help=False,
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128)
    )

    # add the arguments
    parser.add_argument("-g", "--guess", help="Wavefunction guess. (default: %(default)s)", type=str, default="exp(-(x-0.5)**2)")
    parser.add_argument("-h", "--help", help="Print this help message.", action="store_true")
    parser.add_argument("-i", "--iters", help="Maximum number of iterations. (default: %(default)s)", type=int, default=1000)
    parser.add_argument("-n", "--nstate", help="Number of states to consider. (default: %(default)s)", type=int, default=3)
    parser.add_argument("-p", "--points", help="Number of discretization points. (default: %(default)s)", type=int, default=1024)
    parser.add_argument("-r", "--range", help="Range for the calculated wavefunction. (default: %(default)s)", nargs=2, type=float, default=[-16, 16])
    parser.add_argument("-s", "--tstep", help="Time step of the propagation. (default: %(default)s)", type=float, default=0.1)
    parser.add_argument("-t", "--threshold", help="Convergence threshold for the wavefunction. (default: %(default)s)", type=float, default=1e-12)
    parser.add_argument("-v", "--potential", help="Model potential. (default: %(default)s)", type=str, default="0.5*x**2")

    # add the flags
    parser.add_argument("--optimize", help="Enable initial optimization for real time propagation.", action="store_true")
    parser.add_argument("--real", help="Perform the real time propagation.", action="store_true")

    # parse the arguments
    args = parser.parse_args()

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # VARIABLE INITIALIZATION ==========================================================================================================================================================================

    # define the discretization step and state array
    dx, states = (args.range[1] - args.range[0]) / (args.points - 1), list()

    # define x and k space
    x, k = np.linspace(args.range[0], args.range[1], args.points), 2 * np.pi * np.fft.fftfreq(args.points, dx)

    # define initial wavefunction and potential
    psi0, V = eval(args.guess.replace("exp", "np.exp")), eval(args.potential.replace("exp", "np.exp"))

    # define R and K operators for real and imaginary time propagation
    Rr, Kr = np.exp(-0.5j * V * args.tstep), np.exp(-0.5j * k**2 * args.tstep)
    Ri, Ki = np.exp(-0.5 * V * args.tstep), np.exp(-0.5 * k**2 * args.tstep)

    # WAVEFUNCTION PROPAGATION =========================================================================================================================================================================

    # propagate each state
    for i in range(args.nstate):
        # define initial wavefunction
        psi = [0j + psi0 / np.sqrt(np.sum(np.abs(psi0)**2) * dx)]

        # loop over imaginary and real time
        for R, K in list(zip([Ri, Rr], [Ki, Kr])):
            # break if imaginary time propagation is not required
            if args.real and not args.optimize and np.array_equal(R, Ri) and np.array_equal(K, Ki): continue

            # break if real time propagation is not required
            if not args.real and np.array_equal(R, Rr) and np.array_equal(K, Kr): break

            # reset the wavefunction guess for real time propagation
            if np.array_equal(R, Rr) and np.array_equal(K, Kr): psi = [psi[-1]]
            
            # propagate the wavefunction
            for _ in range(args.iters):
                # apply R and K operators and append to wavefunction list
                psi.append(R * np.fft.ifft(K * np.fft.fft(R * psi[-1])));

                # apply Gram-Schmidt orthogonalization
                for j in (j for j in range(i)):
                    psi[-1] -= np.sum(np.conj(states[j][-1]) * psi[-1]) * states[j][-1] * dx

                # normalize wavefunction
                psi[-1] /= np.sqrt(np.sum(np.abs(psi[-1])**2) * dx)

                # break if wavefunction has converged
                if np.sum(np.abs(psi[-1] - psi[-2])**2) < args.threshold or np.abs(energy(psi[-1]) - energy(psi[-2])) < args.threshold: break

        # append wavefunction to list of states and print energy
        states.append(psi); print("E_{}:".format(i), energy(psi[-1]))

    # RESULTS AND PLOTTING =============================================================================================================================================================================

    # define probability density
    D = [[energy(psi) + np.abs(psi)**2 for psi in state] for state in states]

    # create the figure and definte tight layout
    [fig, ax] = plt.subplots(); plt.tight_layout()

    # define minimum and maximum x values for plotting
    xmin = np.min([np.min([np.min(x[np.abs(psi)**2 > 1e-8]) for psi in Si]) for Si in states])
    xmax = np.max([np.max([np.max(x[np.abs(psi)**2 > 1e-8]) for psi in Si]) for Si in states])

    # define maximum y values for plotting
    ymaxreal = max([max([energy(psi) + np.real(psi).max() for psi in state]) for state in states])
    ymaximag = max([max([energy(psi) + np.imag(psi).max() for psi in state]) for state in states])

    # define the minimum y values for plotting
    yminreal = min([min([energy(psi) + np.real(psi).min() for psi in state]) for state in states])
    yminimag = min([min([energy(psi) + np.imag(psi).min() for psi in state]) for state in states])

    # set limits of the plot
    ax.set_xlim(xmin, xmax); ax.set_ylim(np.block([V, yminreal, yminimag]).min(), max([ymaxreal, ymaximag]))

    # plot the potential and initial wavefunctions
    ax.plot(x, V); plots = [[ax.plot(x, np.real(state[0]))[0], ax.plot(x, np.imag(state[0]))[0]] for state in states]

    # animation update function
    def update(j):
        for i in range(len(plots)): plots[i][0].set_ydata(energy(states[i][j if j < len(states[i]) else -1]) + np.real(states[i][j if j < len(states[i]) else -1]))
        for i in range(len(plots)): plots[i][1].set_ydata(energy(states[i][j if j < len(states[i]) else -1]) + np.imag(states[i][j if j < len(states[i]) else -1]))

    # animate the wavefunction
    ani = anm.FuncAnimation(fig, update, frames=np.max([len(state) for state in states]), interval=30); plt.show(); plt.close("all") # type: ignore
