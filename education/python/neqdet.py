#!/usr/bin/env python

import argparse as ap, itertools as it, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, re, scipy as sp, scipy.linalg; np.random.seed(0)

# EXAMPLES
# ./neqdet.py -g "[np.exp(-(r1-1)**2-(r2-1)**2)]" -v "[[0.5*(r1**2+r2**2)]]"
# ./neqdet.py -d 0 -f 0.01 -p 8192 -l 32 -m 2000 -t 10 -v "[[0.01*np.tanh(0.6*r1),0.001*np.exp(-r1**2)],[0.001*np.exp(-r1**2),-0.01*np.tanh(0.6*r1)]]" -g "[0,np.exp(-(r1+10)**2+10j*r1)]" --align --adiabatic

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
    parser.add_argument("-n", "--ntraj", help="Number of Bohmian trajectories. (default: %(default)s)", type=int, default=100)
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

    # calculate the space step, extract the number of dimensions
    dr = 2 * args.limit / (args.points - 1); ndim = max(map(int, re.findall(r"r\[:, (\d)\]", args.potential))) + 1;

    # create the grid in real and fourier space
    r = np.stack(np.meshgrid(*[np.linspace(-args.limit, args.limit, args.points)] * ndim, indexing="ij"), axis=-1).reshape(-1, ndim)
    k = np.stack(np.meshgrid(*[2 * np.pi * np.fft.fftfreq(args.points, dr)      ] * ndim, indexing="ij"), axis=-1).reshape(-1, ndim)

    # define the potential
    V = np.array(eval(args.potential)).transpose(2, 0, 1)

    # REAL AND IMAGINARY QUANTUM DYNAMICS =============================================================================================================================================================

    # create the containers for the observables and wavefunctions
    position, momentum, ekin, epot = [], [], [], []; wfnopt, wfn = [], []

    # iterate over the propagations
    for i in range(args.imaginary if args.imaginary else 1):

        # print the propagation header
        print() if i else None; print("PROPAGATION OF STATE %d " % (i))

        # create the initial wavefunction from the provided guess and normalize it
        psi = np.array(list(map(lambda x: x * np.ones(r.shape[0]), eval(args.guess))), dtype=complex).T; psi /= np.sqrt(dr * np.einsum("ij,ij->", psi.conj(), psi))

        # create the initial points for Bohmian trajectories
        trajs = np.concatenate((r[np.random.choice(psi.shape[0], size=args.ntraj, p=(np.abs(psi)**2 * dr).sum(axis=1))][:, None, :], np.zeros((args.ntraj, args.iterations, ndim))), axis=1)

        # get the full wavefunction shape and clear the containers
        shape = ndim * [args.points] + [psi.shape[1]]; wfn.clear(), position.clear(), momentum.clear(), ekin.clear(), epot.clear()

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
            if (j): psi = np.fft.fftn(psi.reshape(shape), axes=range(ndim)).reshape(psi.shape)

            # propagate in momentum space
            if (j): psi = np.einsum("ijk,ik->ij", K, psi)

            # inverse fourier transform the wavefunction
            if (j): psi = np.fft.ifftn(psi.reshape(shape), axes=range(ndim)).reshape(psi.shape)

            # propagate in real space
            if (j): psi = np.einsum("ijk,ik->ij", R, psi)

            # orthogonalize the wavefunction
            for i in range(len(wfnopt)): psi -= np.sum(wfnopt[i].conj() * psi) * wfnopt[i] * dr

            # normalize the wavefunction
            if args.imaginary: psi /= np.sqrt(dr * np.einsum("ij,ij->", psi.conj(), psi))

            # append the potential energy and the wavefunction
            epot.append(np.einsum("ij,ijk,ik->", psi.conj(), V, psi).real * dr); wfn.append(psi.copy())

            # create a n dimensional copy of the wavefunction, its fourier transform and a container for Bohmian trajectory velocity
            psid, psik, v = psi.reshape(shape), np.fft.fftn(psi.reshape(shape), axes=range(ndim)).reshape(psi.shape), np.zeros(shape[:-1] + [ndim])

            # append the kinetic energy
            ekin.append((psi.conj() * np.fft.ifftn((psik * (0.5 * np.sum(k**2, axis=1) / args.mass)[:, None]).reshape(shape), axes=range(ndim)).reshape(psi.shape)).real.sum() * dr)

            # append the position
            position.append(np.sum(r * np.sum(np.abs(psi)**2, axis=1, keepdims=True), axis=0) * dr)

            # append the momentum
            momentum.append(np.array([np.sum(psi.conj() * np.fft.ifftn((1j * (k[:, dim:dim + 1]) * psik).reshape(shape), axes=range(ndim)).reshape(psi.shape)).imag * dr for dim in range(ndim)]))

            # calculate the velocity of the Bohmian trajectories
            if (j): v[..., :] = np.array([(np.conjugate(psid) * np.gradient(psid, dr, axis=dim)).sum(axis=-1).imag / ((np.abs(psid)**2).sum(axis=-1) + 1e-14) / args.mass for dim in range(ndim)]).T

            # propagate the Bohmian trajectories
            if (j): trajs[:, j, :] = trajs[:, j - 1, :] + sp.interpolate.interpn(points=ndim * [np.unique(r[:, 0])], values=v, xi=trajs[:, j - 1, :]) * args.timestep

            # print the iteration info
            if j % 100 == 0: print("%6d %12.6f %12.6f %12.6f" % (j, ekin[-1], epot[-1], ekin[-1] + epot[-1]))

        # append the optimized wavefunction to the container
        if args.imaginary: wfnopt.append(psi.copy())

    # ADIABATIC TRANSFORM AND SPECTRUM ====================================================================================================================================================================

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

    # scale the wavefunction and add the potential to the wavefunction
    wfn = args.factor * np.array(wfn); wfn += (1 + 1j) * np.einsum("ijj,k->kij", V, np.ones(args.iterations + 1)) if args.align else 0

    # create the stationary subplots and scale the wavefunction
    fig, axs = plt.subplots(2, 3, figsize=(12, 6))

    # plot the energies
    axs[0, 0].plot(np.arange(args.iterations + 1) * args.timestep, ekin, label="Kinetic Energy"  )
    axs[0, 0].plot(np.arange(args.iterations + 1) * args.timestep, epot, label="Potential Energy")

    # plot the population
    axs[0, 1].plot(np.arange(args.iterations + 1) * args.timestep, [np.diag(rho) for rho in density], label=[f"S$_{i}$" for i in range(density.shape[1])])

    # print the acf
    axs[1, 0].plot(np.arange(-args.iterations, args.iterations + 1) * args.timestep, np.array(acf).real, label="Re(ACF)")
    axs[1, 0].plot(np.arange(-args.iterations, args.iterations + 1) * args.timestep, np.array(acf).imag, label="Im(ACF)")

    # plot the spectrum
    axs[1, 1].plot(omega[np.argsort(omega)], spectrum[np.argsort(omega)] / np.max(spectrum))

    # plot the position and momentum of the Bohmian trajectories
    [axs[0, 2].plot(np.arange(args.iterations + 1) * args.timestep, trajs[i, :, 0],                                                 alpha=0.05, color="tab:blue") for i in range(args.ntraj)]
    [axs[1, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.gradient(trajs[i, :, 0], args.timestep, axis=0) * args.mass, alpha=0.05, color="tab:blue") for i in range(args.ntraj)]

    # plot the numerically exact expectation values of position and momentum
    axs[0, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.array(position)[:, 0], color="tab:orange", label="$<\Psi|\hat{r_x}|\Psi>$")
    axs[1, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.array(momentum)[:, 0], color="tab:orange", label="$<\Psi|\hat{p_x}|\Psi>$")

    # plot the mean of the Bohmian trajectories
    axs[0, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.mean(trajs[:, :, 0], axis=0),                                                 "--", color="black", label="$<r_x>$")
    axs[1, 2].plot(np.arange(args.iterations + 1) * args.timestep, np.mean(np.gradient(trajs[:, :, 0], args.timestep, axis=1), axis=0) * args.mass, "--", color="black", label="$<p_x>$")

    # set the labels for the stationary plot
    axs[0, 0].set_xlabel("Time (a.u.)"    ); axs[0, 0].set_ylabel("Energy (a.u.)"           )
    axs[0, 1].set_xlabel("Time (a.u.)"    ); axs[0, 1].set_ylabel("Population"              )
    axs[0, 2].set_xlabel("Time (a.u.)"    ); axs[0, 2].set_ylabel("Position (a.u.)"         )
    axs[1, 0].set_xlabel("Time (a.u.)"    ); axs[1, 0].set_ylabel("Autocorrelation Function")
    axs[1, 1].set_xlabel("Energy (a.u.)"  ); axs[1, 1].set_ylabel("Normalized Intensity"    )
    axs[1, 2].set_xlabel("Time (a.u.)"    ); axs[1, 2].set_ylabel("Momentum (a.u.)"         )

    # set the domain for the spectrum plot, the end will be as last element from the end less than some value
    axs[1, 1].set_xlim(0, omega[np.argsort(omega)][np.where(spectrum[np.argsort(omega)] / np.max(spectrum) > 1e-6)][-1])

    # enable legends
    axs[0, 0].legend(); axs[0, 1].legend(); axs[1, 0].legend(); axs[0, 2].legend(); axs[1, 2].legend()

    # only for 1D
    if ndim == 1:
        
        # creace the wavefunction plot
        wfig, wax = plt.subplots(1, 1)

        # plot the wavefunction
        wfnplot = np.array([[wax.plot(r, wfn[0][:, i].real, label="Re($\Psi$)")[0], wax.plot(r, wfn[0][:, i].imag, label="Im($\Psi$)")[0]] for i in range(psi.shape[1])]).flatten()

        # set the labels for the wavefunction plot and enable legend
        wax.set_xlabel("Position (a.u.)"); wax.set_ylabel("Wavefunction"); wax.legend()

        # extract the wfn min and max
        minwfn, maxwfn = min(np.array(wfn).real.min(), np.array(wfn).imag.min()), max(np.array(wfn).real.max(), np.array(wfn).imag.max())

        # set the limits for the wavefunction animation
        wax.set_ylim(minwfn - 0.1 * (maxwfn - minwfn), maxwfn + 0.1 * (maxwfn - minwfn))

        # define the update function for the wavefunction animation
        update = lambda i: [wfnplot[j].set_ydata(wfn[i][:, j // 2].real if j % 2 == 0 else wfn[i][:, j // 2].imag) for j in range(2 * psi.shape[1])]

        # make the wavefunction plot animation and set the layout
        ani = anm.FuncAnimation(wfig, update, frames=range(len(wfn)), repeat=True, interval=30); wfig.tight_layout()

    # show the plot
    fig.tight_layout(); plt.show()
