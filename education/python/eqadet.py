#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np

def energy(wfn):
    Ek = 0.5 * np.conj(wfn) * np.fft.ifftn(sum([[k, l, m][i]**2 for i in range(dim)]) * np.fft.fftn(wfn))
    Er = np.conj(wfn) * V * wfn; return np.sum(Ek + Er).real * dr

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
    parser.add_argument("-m", "--mass", help="Mass parameter. (default: %(default)s)", type=int, default=1)
    parser.add_argument("-n", "--nstate", help="Number of states to consider. (default: %(default)s)", type=int, default=1)
    parser.add_argument("-p", "--points", help="Number of discretization points. (default: %(default)s)", type=int, default=128)
    parser.add_argument("-r", "--range", help="Range for the calculated wavefunction. (default: %(default)s)", nargs=2, type=float, default=[-16, 16])
    parser.add_argument("-s", "--tstep", help="Time step of the propagation. (default: %(default)s)", type=float, default=0.1)
    parser.add_argument("-t", "--threshold", help="Convergence threshold for the wavefunction. (default: %(default)s)", type=float, default=1e-12)
    parser.add_argument("-v", "--potential", help="Potential where the optimization happens and, if -e not provided, also real time propagation. (default: %(default)s)", type=str, default="0.5*x**2")
    parser.add_argument("-e", "--excpotential", help="Excited state potential where the real time propagation happens. (default: %(default)s)", type=str)

    # add the flags
    parser.add_argument("--optimize", help="Enable initial optimization for real time propagation.", action="store_true")
    parser.add_argument("--static", help="Make the optimization animation show only the last step.", action="store_true")
    parser.add_argument("--real", help="Perform the real time propagation.", action="store_true")

    # parse the arguments
    args = parser.parse_args()

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # VARIABLE INITIALIZATION ==========================================================================================================================================================================

    # define the discretization step, state array and iteration break flag
    dx, states = (args.range[1] - args.range[0]) / (args.points - 1), [[] for _ in range(args.nstate)]; dy, dz = dx, dx;

    # define x and k spaces
    x, k = np.linspace(args.range[0], args.range[1], args.points), 2 * np.pi * np.fft.fftfreq(args.points, dx); y, z = x, x; l, m = k, k;

    # define the number of dimensions and coordinate step
    dim = 3 if "z" in args.potential else (2 if "y" in args.potential else 1); dr = np.prod([[dx, dy, dz][i] for i in range(dim)])

    # redefine the variables for more dimensions
    if dim == 3: x, y, z = np.meshgrid(x, y, z); k, l, m = np.meshgrid(k, l, m)
    if dim == 2: x, y = np.meshgrid(x, y); k, l = np.meshgrid(k, l)

    # define the time and frequency axis
    t, f = np.linspace(0, args.iters * args.tstep, args.iters + 1), 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(args.iters + 1, args.tstep))

    # define initial wavefunction and potential
    psi0, V = eval(args.guess.replace("exp", "np.exp")), eval(args.potential.replace("exp", "np.exp"))

    # WAVEFUNCTION PROPAGATION =========================================================================================================================================================================

    # define R and K operators for imaginary time propagation
    R, K = np.exp(-0.5 * V * args.tstep), np.exp(-0.5 * sum([[k, l, m][i]**2 for i in range(dim)]) * args.tstep / args.mass)

    # propagate each state in imaginary time if requested
    for i in (i for i in range(args.nstate) if not args.real or (args.real and args.optimize)):
        # define initial wavefunction
        psi = [0j + psi0 / np.sqrt(np.sum(np.abs(psi0)**2) * dr)]

        # propagate the wavefunction
        for _ in range(args.iters):
            # apply R and K operators and append to wavefunction list
            psi.append(R * np.fft.ifftn(K * np.fft.fftn(R * psi[-1])));

            # apply Gram-Schmidt orthogonalization
            for j in (j for j in range(i)):
                psi[-1] -= np.sum(np.conj(states[j][-1]) * psi[-1]) * states[j][-1] * dr

            # normalize wavefunction
            psi[-1] /= np.sqrt(np.sum(np.abs(psi[-1])**2) * dr)

        # append wavefunction to list of states and print energy
        states[i] = psi; print("E_{}:".format(i), energy(psi[-1]))

    # change the potential to the excited one if provided
    if args.excpotential: V = eval(args.excpotential.replace("exp", "np.exp"))

    # define R and K operators for real time propagation
    R, K = np.exp(-0.5j * V * args.tstep), np.exp(-0.5j * k**2 * args.tstep / args.mass)

    # propagate each state in real time if requested
    for i in (i for i in range(args.nstate) if args.real):
        # define initial wavefunction
        psi = [states[i][-1]] if args.optimize else [0j + psi0 / np.sqrt(np.sum(np.abs(psi0)**2) * dr)]

        # propagate the wavefunction
        for _ in range(args.iters): psi.append(R * np.fft.ifft(K * np.fft.fft(R * psi[-1])));

        # append wavefunction to list of states and print energy
        states[i] = psi; print("E_{}:".format(i), energy(psi[-1]))

    # calculate the autocorrelation function of a ground state and it Fourier transform
    if args.real: G = np.array([np.sum((states[0][0]) * np.conj(psi)) * dr for psi in states[0]]); F = np.fft.fftshift(np.fft.fftn(G))

    # RESULTS AND PLOTTING =============================================================================================================================================================================

    # start frame
    sf = -1 if args.static else 0

    if dim == 1:
        # create the figure and definte tight layout
        [fig, ax] = plt.subplots(1, 2 if args.real else 1, figsize=(12, 5) if args.real else (6, 5)); ax = ax if args.real else [ax]; plt.tight_layout()

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
        ax[0].set_xlim(xmin, xmax); ax[0].set_ylim(np.block([V, yminreal, yminimag]).min(), max([ymaxreal, ymaximag]))

        # plot the potential and initial wavefunctions
        ax[0].plot(x, V); plots = [[ax[0].plot(x, energy(state[sf]) + np.real(state[sf]))[0], ax[0].plot(x, energy(state[sf]) + np.imag(state[sf]))[0]] for state in states]

        # animation update function
        def update(j):
            for i in range(len(plots)): plots[i][0].set_ydata(energy(states[i][j if j < len(states[i]) else -1]) + np.real(states[i][j if j < len(states[i]) else -1]))
            for i in range(len(plots)): plots[i][1].set_ydata(energy(states[i][j if j < len(states[i]) else -1]) + np.imag(states[i][j if j < len(states[i]) else -1]))

        # plot the spectrum
        if args.real: ax[1].plot(f, np.abs(F)) # type: ignore

    if dim == 2:
        # create the figure
        [fig, ax] = plt.subplots((args.nstate - 1) // 5 + 1, min([5, args.nstate]), figsize=(3.2 * min([5, args.nstate]), ((args.nstate - 1) // 5 + 1) * 3))

        # define spacing and modify axes if flat
        plt.tight_layout(); plt.subplots_adjust(wspace=0.15); ax = ax if args.nstate > 1 else np.array([ax])

        # define the plots
        plots = [ax.pcolormesh(x, y, np.abs(states[i][sf])) for i, ax in enumerate(ax.flat) if i < args.nstate]

        # set axes limits
        [[ax.set_xlim(-6, 6), ax.set_ylim(-6, 6)] for _, ax in enumerate(ax.flat)]

        # animation update function
        def update(j):
            for i in range(len(plots)): plots[i].set_array(np.abs(states[i][j if j < len(states[i]) else -1]))
        
    # animate the wavefunction
    if not args.static or args.real: ani = anm.FuncAnimation(fig, update, frames=np.max([len(state) for state in states]), interval=30) # type: ignore

    # set the window name
    fig.canvas.manager.set_window_title("Exact Quantum Adiabatic Dynamics Educational Toolkit") # type: ignore

    # shot the plots
    plt.show(); plt.close("all")
