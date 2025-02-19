#!/usr/bin/env python

import argparse as ap, itertools as it, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, re, scipy as sp, scipy.linalg

# EXAMPLES
# ./neqdet.py -n 2 -w "[np.exp(-(r1-1)**2-(r2-1)**2)]" -v "[[0.5*(r1**2+r2**2)]]"
# ./neqdet.py -d 0 -f 0.01 -p 2048 -l 32 -m 2000 -t 10 -v "[[0.01*np.tanh(0.6*r1),0.001*np.exp(-r1**2)],[0.001*np.exp(-r1**2),-0.01*np.tanh(0.6*r1)]]" -g "[0,np.exp(-(r1+10)**2+10j*r1)]" --align

if __name__ == "__main__":
    # SECTION FOR PARSING COMMAND LINE ARGUMENTS =======================================================================================================================================================

    # create the parser
    parser = ap.ArgumentParser(
        prog="resmet.py", description="Numerically Exact Quantum Dynamics Educational Toolkit",
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    # add the arguments
    parser.add_argument("-a", "--align", help="Align the wavefunction plot to the potential.", action=ap.BooleanOptionalAction)
    parser.add_argument("-c", "--imaginary", help="Perform imaginary-time propagation and find specified number of orthogonal states.", type=int, default=0)
    parser.add_argument("-d", "--damp", help="Gaussian damping parameter for autocorrelation function. (default: %(default)s)", type=float, default=0.003)
    parser.add_argument("-f", "--factor", help="Factor to scale the wavefunction before plotting. (default: %(default)s)", type=float, default=10)
    parser.add_argument("-g", "--guess", help="Initial guess for the wavefunction. (default: %(default)s)", type=str, default="[np.exp(-(r1-1)**2)]")
    parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)
    parser.add_argument("-i", "--iterations", help="Number of iterations. (default: %(default)s)", type=int, default=500)
    parser.add_argument("-l", "--limit", help="Distance limit of the wavefunction grid. (default: %(default)s)", type=float, default=8)
    parser.add_argument("-m", "--mass", help="Mass of the particle. (default: %(default)s)", type=float, default=1)
    parser.add_argument("-p", "--points", help="Number of points in the wavefunction grid. (default: %(default)s)", type=int, default=128)
    parser.add_argument("-t", "--timestep", help="Time step of the simulation. (default: %(default)s)", type=float, default=0.1)
    parser.add_argument("-u", "--adiabatic", help="Transform the results to the adiabatic basis. (default: %(default)s)", action=ap.BooleanOptionalAction)
    parser.add_argument("-v", "--potential", help="Potential matrix. (default: %(default)s)", type=str, default="[[0.5*r1**2]]")

    # parse the arguments
    args = parser.parse_args()

    # replace the variables in the string expressions
    for i in range(1, 10): args.guess     =     args.guess.replace(f"r{i}", f"r[:, {i - 1}]")
    for i in range(1, 10): args.potential = args.potential.replace(f"r{i}", f"r[:, {i - 1}]")

    # print the help message if the flag is set
    if args.help: parser.print_help(); exit()

    # PREPARE VARIABLES FOR THE DYNAMICS ==============================================================================================================================================================

    # calculate the space step, extract the number of dimensions and create the containers for the wavefunctions
    dr = 2 * args.limit / (args.points - 1); ndim = max(map(int, re.findall(r"r\[:, (\d)\]", args.potential))) + 1; wfnopt, wfn = [], []

    # create the grid in real and fourier space
    r = np.stack(np.meshgrid(*[np.linspace(-args.limit, args.limit, args.points)] * ndim, indexing="ij"), axis=-1).reshape(-1, ndim)
    k = np.stack(np.meshgrid(*[2 * np.pi * np.fft.fftfreq(args.points, dr)      ] * ndim, indexing="ij"), axis=-1).reshape(-1, ndim)

    # define the potential
    V = np.array(eval(args.potential)).transpose(2, 0, 1)

    # PERFORM THE PROPAGATIONS =========================================================================================================================================================================

    # iterate over the propagations
    for i in range(args.imaginary if args.imaginary else 1):

        # print the propagation header and clear the wfn container
        print("%sPROPAGATION OF STATE %d " % ("\n " if i else "", i)); wfn.clear()

        # create the initial wavefunction from the provided guess
        psi = np.array(list(map(lambda x: x*np.ones(r.shape[0]), eval(args.guess))), dtype=complex).T; psi0 = psi.copy() / np.sqrt(np.sum(psi.conj() * psi) * dr)

        # calculate the propagators for each point in the grid
        K = np.array([sp.linalg.expm(-0.5 * (1 if args.imaginary else 1j) * args.timestep * np.sum(k[i, :] ** 2) / args.mass * np.eye(psi.shape[1])) for i in range(r.shape[0])])
        R = np.array([sp.linalg.expm(-0.5 * (1 if args.imaginary else 1j) * args.timestep                                    * V[i]                ) for i in range(r.shape[0])])

        # print the propagation header
        print("%6s %12s %12s %12s" % ("ITER", "EKIN", "EPOT", "ETOT", ))

        # propagate the wavefunction
        for j in range(args.iterations + 1):

            # propagate in real space 
            if (j): psi = np.einsum("ijk,ik->ij", R, psi)

            # fourier transform the wavefunction
            if (j): psi = np.fft.fftn(psi.reshape(ndim * [args.points] + [psi.shape[1]]), axes=range(ndim)).reshape(psi.shape)

            # propagate in momentum space
            if (j): psi = np.einsum("ijk,ik->ij", K, psi)

            # inverse fourier transform the wavefunction
            if (j): psi = np.fft.ifftn(psi.reshape(ndim * [args.points] + [psi.shape[1]]), axes=range(ndim)).reshape(psi.shape)

            # propagate in real space
            if (j): psi = np.einsum("ijk,ik->ij", R, psi)

            # orthogonalize the wavefunction
            for i in range(len(wfnopt)): psi -= np.sum(wfnopt[i].conj() * psi) * wfnopt[i] * dr

            # normalize the wavefunction
            if args.imaginary or j == 0: psi /= np.sqrt(dr * np.einsum("ij,ij->", psi.conj(), psi))

            # calculate the potential energy and append the wavefunction
            Epot = np.einsum("ij,ijk,ik->", psi.conj(), V, psi).real * dr; wfn.append(psi.copy())

            # create a copy of the wavefunction and store there the fourier transform of the wavefunction
            psik = np.fft.fftn(psi.reshape(ndim * [args.points] + [psi.shape[1]]), axes=range(ndim)).reshape(psi.shape)

            # multiply each state of the transformed wavefunction by k
            for i in range(psi.shape[1]): psik[:, i] *= 0.5 * np.sum(k**2, axis=1) / args.mass

            # calculate the kinetic energy
            Ekin = np.sum(psi.conj() * np.fft.ifftn(psik.reshape(ndim * [args.points] + [psi.shape[1]]), axes=range(ndim)).reshape(psi.shape)).real * dr

            # print the iteration info
            if j % 100 == 0: print("%6d %12.6f %12.6f %12.6f" % (j, Ekin, Epot, Ekin + Epot))

        # append the optimized wavefunction to the container
        if args.imaginary: wfnopt.append(psi.copy())

    # calculate the adiabatic eigenstates
    U = [np.linalg.eigh(V[i])[1] for i in range(V.shape[0])]

    # adiabatize the potential and wavefunctions
    V   = np.einsum("ijk,ikl,ilm->ijm", U, V,   U) if args.adiabatic else np.array(V  )
    wfn = np.einsum("ikl,jik->jil",     U, wfn   ) if args.adiabatic else np.array(wfn)

    # calculate the density matrices and acf
    density = np.einsum("jia,jib->jab", wfn,           wfn.conj()).real * dr
    acf     = np.einsum("ij,tij->t",    wfn[0].conj(), wfn       )      * dr

    # symmetrize the acf and apply the damping function
    acf = np.concatenate((np.flip(acf)[:-1], np.array(acf).conj())) * np.exp(-args.damp * (np.arange(-args.iterations, args.iterations + 1) * args.timestep)**2)

    # calculate the spectrum of the zero-padded acf and the corresponding energies
    spectrum = np.abs(np.fft.fft(np.pad(acf, 2 * [10 * len(acf)], mode="constant")))**2; omega = 2 * np.pi * np.fft.fftfreq(len(spectrum), args.timestep)

    # PRINT AND PLOT THE RESULTS ==========================================================================================================================================================================

    # print the final population
    print("\nFINAL POPULATION: %s" % np.diag(density[-1]))

    # create the subplots and scale the wavefunction
    fig, axs = plt.subplots(2, 2, figsize=(8, 6)); wfn = args.factor * np.array(wfn)

    # add the potential to the wavefunction
    if args.align: wfn += (1 + 1j) * np.einsum("ijj,k->kij", V, np.ones(args.iterations + 1))

    # plot the wavefunction
    if ndim == 1: wfnplot = np.array([[axs[0, 0].plot(r, wfn[0][:, i].real)[0], axs[0, 0].plot(r, wfn[0][:, i].imag)[0]] for i in range(psi.shape[1])]).flatten()

    # plot the population
    axs[0, 1].plot(np.arange(args.iterations + 1) * args.timestep, [np.diag(rho) for rho in density])

    # print the acf
    axs[1, 0].plot(np.arange(-args.iterations, args.iterations + 1) * args.timestep, np.array(acf).real)
    axs[1, 0].plot(np.arange(-args.iterations, args.iterations + 1) * args.timestep, np.array(acf).imag)

    # plot the spectrum
    axs[1, 1].plot(omega[np.argsort(omega)], spectrum[np.argsort(omega)] / np.max(spectrum))

    # set the labels
    axs[0, 0].set_xlabel("Position (a.u.)"); axs[0, 0].set_ylabel("Wavefunction"            )
    axs[0, 1].set_xlabel("Time (a.u.)"    ); axs[0, 1].set_ylabel("Population"              )
    axs[1, 0].set_xlabel("Time (a.u.)"    ); axs[1, 0].set_ylabel("Autocorrelation Function")
    axs[1, 1].set_xlabel("Energy (a.u.)"  ); axs[1, 1].set_ylabel("Normalized Intensity"    )

    # extract the wfn min and max
    minwfn, maxwfn = min(np.array(wfn).real.min(), np.array(wfn).imag.min()), max(np.array(wfn).real.max(), np.array(wfn).imag.max())

    # set the domain for the spectrum plot, the end will be as last element from the end less than some value
    axs[1, 1].set_xlim(0, omega[np.argsort(omega)][np.where(spectrum[np.argsort(omega)] / np.max(spectrum) > 1e-6)][-1])

    # set the limits for the wavefunction animation
    axs[0, 0].set_ylim(minwfn - 0.1 * (maxwfn - minwfn), maxwfn + 0.1 * (maxwfn - minwfn))

    # define the update function for the wavefunction animation
    update = lambda i: [wfnplot[j].set_ydata(wfn[i][:, j // 2].real if j % 2 == 0 else wfn[i][:, j // 2].imag) for j in range(2 * psi.shape[1])]

    # make the wavefunction plot animation
    if ndim == 1: ani = anm.FuncAnimation(fig, update, frames=range(len(wfn)), repeat=True, interval=30)

    # show the plot
    plt.tight_layout(); plt.show()
