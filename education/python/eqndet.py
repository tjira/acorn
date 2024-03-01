#!/usr/bin/env python

import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np

def potexp(V):
    D = 4 * np.abs(V[1, 0])**2 + (V[0, 0] - V[1, 1])**2
    R = np.exp(-0.25j * (V[0, 0] + V[1, 1]) * args.tstep) * (np.array([[np.cos(0.25 * np.sqrt(D) * args.tstep), 0 * V[0, 0]], [0 * V[0, 0], np.cos(0.25 * np.sqrt(D) * args.tstep)]])
      + 1.0j * np.sin(0.25 * np.sqrt(D) * args.tstep) / np.sqrt(D) * np.array([[V[1, 1] - V[0, 0], -2 * V[0, 1]], [-2 * V[1, 0], V[0, 0] - V[1, 1]]]))
    return R

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =======================================================================================================================================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="eqdet.py", description="Exact Quantum Adiabatic Dynamics Educational Toolkit", add_help=False,
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128)
    )

    # add the arguments
    parser.add_argument("-g", "--guess", help="Wavefunction guess. (default: %(default)s)", type=str, default="exp(-(x+15)**2)")
    parser.add_argument("-h", "--help", help="Print this help message.", action="store_true")
    parser.add_argument("-i", "--iters", help="Maximum number of iterations. (default: %(default)s)", type=int, default=2048)
    parser.add_argument("-m", "--mass", help="Mass parameter. (default: %(default)s)", type=int, default=2048)
    parser.add_argument("-p", "--points", help="Number of discretization points. (default: %(default)s)", type=int, default=2048)
    parser.add_argument("-r", "--range", help="Range for the calculated wavefunction. (default: %(default)s)", nargs=2, type=float, default=[-32, 32])
    parser.add_argument("-s", "--tstep", help="Time step of the propagation. (default: %(default)s)", type=float, default=32)

    # parse the arguments
    args = parser.parse_args()

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # VARIABLE INITIALIZATION ==========================================================================================================================================================================

    # define the discretization step
    dx = (args.range[1] - args.range[0]) / (args.points - 1)

    # define x and k space
    x, k = np.linspace(args.range[0], args.range[1], args.points), 2 * np.pi * np.fft.fftfreq(args.points, dx)

    # define the potential
    V = np.array([
        [+0.01 * np.tanh(0.6 * x), 0.001 * np.exp(-x**2)],
        [0.001 * np.exp(-x**2), -0.01 * np.tanh(0.6 * x)]
    ])

    # define initial wavefunction
    psi = [np.array([eval(args.guess.replace("exp", "np.exp")) * np.exp(1j * np.sqrt(0.06 * args.mass) * x), 0 * x])]; psi[0][0] /= np.sqrt(np.vdot(psi, psi) * dx)

    # WAVEFUNCTION PROPAGATION =========================================================================================================================================================================

    # define R and K operators for real time propagation
    K, R = np.array([[np.exp(-0.5j * k**2 * args.tstep / args.mass), 0 * k], [0 * k, np.exp(-0.5j * k**2 * args.tstep / args.mass)]]), potexp(V)

    # propagate the wavefunction
    for _ in range(args.iters):
        newpsi1 = np.fft.fft(R[0, 0] * psi[-1][0] + R[0, 1] * psi[-1][1])
        newpsi2 = np.fft.fft(R[1, 0] * psi[-1][0] + R[1, 1] * psi[-1][1])
        newpsi1 = np.fft.ifft(K[0, 0] * newpsi1 + K[0, 1] * newpsi2)
        newpsi2 = np.fft.ifft(K[1, 0] * newpsi1 + K[1, 1] * newpsi2)
        newpsi1 = R[0, 0] * newpsi1 + R[0, 1] * newpsi2
        newpsi2 = R[1, 0] * newpsi1 + R[1, 1] * newpsi2
        psi.append(np.array([newpsi1, newpsi2]));

    # RESULTS AND PLOTTING =============================================================================================================================================================================

    # scale the potential
    V[0, 0] *= max([np.real(psi).max(), np.imag(psi).max()]) / V[0, 0].max()
    V[1, 1] *= max([np.real(psi).max(), np.imag(psi).max()]) / V[1, 1].max()

    # create the figure and define tight layout
    [fig, ax] = plt.subplots(); plt.tight_layout()

    # set the window name
    fig.canvas.manager.set_window_title("Exact Quantum Nonadiabatic Dynamics Educational Toolkit") # type: ignore

    # define maximum y values for plotting
    ymaxreal = np.real(psi).max() + V[1, 1].max()
    ymaximag = np.imag(psi).max() + V[1, 1].max()

    # define minimum y values for plotting
    yminreal = np.real(psi).min() + V[0, 0].min()
    yminimag = np.imag(psi).min() + V[0, 0].min()

    # set limits of the plot
    ax.set_ylim(min([yminreal, yminimag]), max([ymaxreal, ymaximag]))

    # plot the potential
    ax.plot(x, V[0, 0]); ax.plot(x, V[1, 1])

    # plot the initial wavefunctions
    plots = [
        ax.plot(x, np.real(psi[0][0]) + V[0, 0].min())[0],
        ax.plot(x, np.imag(psi[0][0]) + V[0, 0].min())[0],
        ax.plot(x, np.real(psi[0][1]) + V[1, 1].max())[0],
        ax.plot(x, np.imag(psi[0][1]) + V[1, 1].max())[0]
    ]

    # animation update function
    def update(i):
        plots[0].set_ydata(np.real(psi[i][0]) + V[0, 0].min())
        plots[1].set_ydata(np.imag(psi[i][0]) + V[0, 0].min())
        plots[2].set_ydata(np.real(psi[i][1]) + V[1, 1].max())
        plots[3].set_ydata(np.imag(psi[i][1]) + V[1, 1].max())

    # animate the wavefunction
    ani = anm.FuncAnimation(fig, update, frames=len(psi), interval=30) # type: ignore

    # shot the plots
    plt.show(); plt.close("all")
