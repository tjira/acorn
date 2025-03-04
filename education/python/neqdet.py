#!/usr/bin/env python

import argparse as ap, itertools as it, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, scipy as sp, scipy.linalg

np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%20.14f" % x)); np.random.seed(0)

# EXAMPLES
# ./neqdet.py -g "[np.exp(-(r1-1)**2-(r2-1)**2)]" -v "[[0.5*(r1**2+r2**2)]]"
# ./neqdet.py -d 0 -f 0.01 -p 8192 -l 32 -m 2000 -t 10 -v "[[0.01*np.tanh(0.6*r1),0.001*np.exp(-r1**2)],[0.001*np.exp(-r1**2),-0.01*np.tanh(0.6*r1)]]" -g "[0,np.exp(-(r1+10)**2+10j*r1)]" --align --adiabatic

# SECTION FOR PARSING COMMAND LINE ARGUMENTS =====================================================================================

# create the parser
parser = ap.ArgumentParser(
    prog="resmet.py", description="Numerically Exact Quantum Dynamics Educational Toolkit",
    formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
    add_help=False, allow_abbrev=False
)

# add the arguments
parser.add_argument("-a", "--align", help="Align the wavefunction plot to the potential.", action=ap.BooleanOptionalAction)
parser.add_argument("-c", "--imaginary", help="Perform imaginary-time propagation for n states.", type=int, default=0)
parser.add_argument("-d", "--damp", help="Gaussian damping parameter. (default: %(default)s)", type=float, default=0.003)
parser.add_argument("-f", "--factor", help="Factor to scale the wavefunction. (default: %(default)s)", type=float, default=10)
parser.add_argument("-g", "--guess", help="Initial guess. (default: %(default)s)", type=str, default="[np.exp(-(r1-1)**2)]")
parser.add_argument("-h", "--help", help="Print this help message.", action=ap.BooleanOptionalAction)
parser.add_argument("-i", "--iterations", help="Number of iterations. (default: %(default)s)", type=int, default=500)
parser.add_argument("-l", "--limit", help="Grid limits. (default: %(default)s)", type=float, default=8)
parser.add_argument("-m", "--mass", help="Mass of the particle. (default: %(default)s)", type=float, default=1)
parser.add_argument("-n", "--ntraj", help="Number of Bohmian trajectories. (default: %(default)s)", type=int, default=100)
parser.add_argument("-p", "--points", help="Number of grid points. (default: %(default)s)", type=int, default=128)
parser.add_argument("-t", "--timestep", help="Time step of the simulation. (default: %(default)s)", type=float, default=0.1)
parser.add_argument("-u", "--adiabatic", help="Adiabatic transform. (default: %(default)s)", action=ap.BooleanOptionalAction)
parser.add_argument("-v", "--potential", help="Potential matrix. (default: %(default)s)", type=str, default="[[0.5*r1**2]]")

# parse the arguments
args = parser.parse_args()

# replace the variables in the string expressions
for i in range(1, 10): args.guess     =     args.guess.replace(f"r{i}", f"r[:, {i - 1}]")
for i in range(1, 10): args.potential = args.potential.replace(f"r{i}", f"r[:, {i - 1}]")

# print the help message if the flag is set
if args.help: parser.print_help(); exit()

# PREPARE VARIABLES FOR THE DYNAMICS =============================================================================================

# extract the number of dimensions without using regex
ndim = max([int(e[0]) for e in args.potential.split("r[:, ")[1:]]) + 1

# define the function to evaluate potential
potf = lambda: np.array(eval(args.potential)).transpose(2, 0, 1);

# define the function to evaluate the guess wavefunction on a specified grid
psif = lambda: np.array(list(map(lambda x: x * np.ones(args.points**ndim), eval(args.guess))), dtype=complex).T

# REAL AND IMAGINARY QUANTUM DYNAMICS ============================================================================================

# calculate the space step and the grid limits in each dimension
dr, grid = 2 * args.limit / (args.points - 1), [np.linspace(-args.limit, args.limit, args.points)] * ndim

# create the grid in real and fourier space
r = np.stack(np.meshgrid(*[np.linspace(-args.limit, args.limit, args.points)] * ndim, indexing="ij"), axis=-1).reshape(-1, ndim)
k = np.stack(np.meshgrid(*[2 * np.pi * np.fft.fftfreq(args.points, dr)      ] * ndim, indexing="ij"), axis=-1).reshape(-1, ndim)

# create the potential, wavefunction and extract wavefunction dimensions
V, psi = potf(), psif(); shape = ndim * [args.points] + [psi.shape[1]]

# normalize the guess wavefunction
psi /= np.sqrt(dr * np.einsum("ij,ij->", psi.conj(), psi))

# create the containers for the observables and wavefunctions
position, momentum, ekin, epot = [], [], [], []; wfnopt, wfn = [], [psi]

# create the function to apply an operator in the momentum space to the n dimensional wavefunction
fapply = lambda op, psik: np.fft.ifftn((op * psik).reshape(shape), axes=range(ndim)).reshape(psi.shape)

# iterate over the propagations
for i in range(args.imaginary if args.imaginary else 1):

    # print the propagation header
    print() if i else None; print("PROPAGATION OF STATE %d " % (i))

    # get the random coordinate indices
    rind = np.random.choice(psi.shape[0], size=args.ntraj, p=(np.abs(psi)**2 * dr).sum(axis=1))

    # create the initial points for Bohmian trajectories
    trajs = np.concatenate((r[rind ][:, None, :], np.zeros((args.ntraj, args.iterations, ndim))), axis=1)

    # initialize the initial wavefunction from the wfn container and clear all containers
    psi = wfn[0]; wfn.clear(), position.clear(), momentum.clear(), ekin.clear(), epot.clear()

    # get the propagator exponent unit
    unit = -0.5 * (1 if args.imaginary else 1j) * args.timestep

    # calculate the propagators for each point in the grid
    K = np.array([sp.linalg.expm(unit * np.sum(k[i, :]**2) / args.mass * np.eye(psi.shape[1])) for i in range(r.shape[0])])
    R = np.array([sp.linalg.expm(unit *                                  V[i]                ) for i in range(r.shape[0])])

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

        # create a n dimensional copy of the wavefunction and its fourier transform
        psid, psik = psi.reshape(shape), np.fft.fftn(psi.reshape(shape), axes=range(ndim)).reshape(psi.shape)

        # append the kinetic energy
        ekin.append((psi.conj() * fapply((0.5 * np.sum(k**2, axis=1) / args.mass)[:, None], psik)).real.sum() * dr)

        # append the position
        position.append(np.sum(r * np.sum(np.abs(psi)**2, axis=1, keepdims=True), axis=0) * dr)

        # append the momentum
        momentum.append([(psi.conj() * fapply(1j * k[:, l:l + 1], psik)).imag.sum() * dr for l in range(ndim)])

        # Bohmian velocity container
        v = np.zeros(shape[:-1] + [ndim])

        # calculate probability current in each dimension
        if (j): currents = [(np.conjugate(psid) * np.gradient(psid, dr, axis=l)).sum(axis=-1).imag for l in range(ndim)]

        # calculate the velocity of the Bohmian trajectories
        if (j): v[..., :] = np.array([currents[l] / ((np.abs(psid)**2).sum(axis=-1) + 1e-14) / args.mass for l in range(ndim)]).T

        # propagate the Bohmian trajectories
        if (j): trajs[:, j] = (lambda r: r + sp.interpolate.interpn(points=grid, values=v, xi=r) * args.timestep)(trajs[:, j - 1])

        # print the iteration info
        if j % 100 == 0: print("%6d %12.6f %12.6f %12.6f" % (j, ekin[-1], epot[-1], ekin[-1] + epot[-1]))

    # append the optimized wavefunction to the container
    if args.imaginary: wfnopt.append(psi.copy())

# ADIABATIC TRANSFORM AND SPECTRUM ===============================================================================================

# calculate the adiabatic eigenstates
U = [np.linalg.eigh(V[i])[1] for i in range(V.shape[0])]

# adiabatize the potential and wavefunctions
V   = np.einsum("ijk,ikl,ilm->ijm", U, V,   U) if args.adiabatic else np.array(V  )
wfn = np.einsum("ikl,jik->jil",     U, wfn   ) if args.adiabatic else np.array(wfn)

# calculate the density matrices and acf
density = np.einsum("jia,jib->jab", wfn,           wfn.conj()).real * dr
acf     = np.einsum("ij,tij->t",    wfn[0].conj(), wfn       )      * dr

# symmetrize the acf
acf = np.concatenate((np.flip(acf)[:-1], np.array(acf).conj()))

# apply the damping to the acf
acf *= np.exp(-args.damp * (np.arange(-args.iterations, args.iterations + 1) * args.timestep)**2)

# create the padded acf
acfpad = np.pad(acf, 2 * [10 * len(acf)], mode="constant")

# calculate the spectrum of the zero-padded acf and the corresponding energies
spectrum = np.abs(np.fft.fft(acfpad))**2; omega = 2 * np.pi * np.fft.fftfreq(len(spectrum), args.timestep)

# PRINT AND PLOT THE RESULTS =====================================================================================================

# print the final population
print(); print("FINAL POPULATION: %s" % np.diag(density[-1]))

# scale the wavefunction and add the potential to the wavefunction
wfn = args.factor * np.array(wfn); wfn += (1 + 1j) * np.einsum("ijj,k->kij", V, np.ones(args.iterations + 1)) if args.align else 0

# extract the momenta from the Bohmian trajectories
ptrajs = np.gradient(trajs, args.timestep, axis=1) * args.mass

# create the stationary subplots
fig, axs = plt.subplots(2, 3, figsize=(12, 6)); time = np.arange(args.iterations + 1) * args.timestep

# plot the energies
axs[0, 0].plot(time, ekin, label="Kinetic Energy"  )
axs[0, 0].plot(time, epot, label="Potential Energy")

# plot the population
axs[0, 1].plot(time, [np.diag(rho) for rho in density], label=[f"S$_{i}$" for i in range(density.shape[1])])

# print the acf
axs[1, 0].plot(np.arange(-args.iterations, args.iterations + 1) * args.timestep, np.array(acf).real, label="Re(ACF)")
axs[1, 0].plot(np.arange(-args.iterations, args.iterations + 1) * args.timestep, np.array(acf).imag, label="Im(ACF)")

# plot the spectrum
axs[1, 1].plot(omega[np.argsort(omega)], spectrum[np.argsort(omega)] / np.max(spectrum))

# plot the position and momentum of the Bohmian trajectories
[axs[0, 2].plot(time,  trajs[i, :, 0], alpha=0.05, color="tab:blue") for i in range(args.ntraj)]
[axs[1, 2].plot(time, ptrajs[i, :, 0], alpha=0.05, color="tab:blue") for i in range(args.ntraj)]

# plot the numerically exact expectation values of position and momentum
axs[0, 2].plot(time, np.array(position)[:, 0], color="tab:orange", label="$<\Psi|\hat{r_x}|\Psi>$")
axs[1, 2].plot(time, np.array(momentum)[:, 0], color="tab:orange", label="$<\Psi|\hat{p_x}|\Psi>$")

# plot the mean of the Bohmian trajectories
axs[0, 2].plot(time, np.mean( trajs[:, :, 0], axis=0), "--", color="black", label="$<r_x>$")
axs[1, 2].plot(time, np.mean(ptrajs[:, :, 0], axis=0), "--", color="black", label="$<p_x>$")

# set the labels for the stationary plot
axs[0, 0].set_xlabel("Time (a.u.)"  ); axs[0, 0].set_ylabel("Energy (a.u.)"           )
axs[0, 1].set_xlabel("Time (a.u.)"  ); axs[0, 1].set_ylabel("Population"              )
axs[0, 2].set_xlabel("Time (a.u.)"  ); axs[0, 2].set_ylabel("Position (a.u.)"         )
axs[1, 0].set_xlabel("Time (a.u.)"  ); axs[1, 0].set_ylabel("Autocorrelation Function")
axs[1, 1].set_xlabel("Energy (a.u.)"); axs[1, 1].set_ylabel("Normalized Intensity"    )
axs[1, 2].set_xlabel("Time (a.u.)"  ); axs[1, 2].set_ylabel("Momentum (a.u.)"         )

# set the domain for the spectrum plot, the end will be as last element from the end less than some value
axs[1, 1].set_xlim(0, omega[np.argsort(omega)][np.where(spectrum[np.argsort(omega)] / np.max(spectrum) > 1e-6)][-1])

# enable legends
axs[0, 0].legend(); axs[0, 1].legend(); axs[1, 0].legend(); axs[0, 2].legend(); axs[1, 2].legend()

# only for 1D
if ndim == 1:
    
    # creace the wavefunction plot
    wfig, wax = plt.subplots(1, 1)

    # plot the wavefunction
    wfnplot = np.array([
        [wax.plot(r, wfn[0][:, i].real, label="Re($\Psi$)")[0], wax.plot(r, wfn[0][:, i].imag, label="Im($\Psi$)")[0]]
    for i in range(psi.shape[1])]).flatten()

    # set the labels for the wavefunction plot and enable legend
    wax.set_xlabel("Position (a.u.)"); wax.set_ylabel("Wavefunction"); wax.legend()

    # extract the wfn min and max
    minwfn, maxwfn = min(np.array(wfn).real.min(), np.array(wfn).imag.min()), max(np.array(wfn).real.max(), np.array(wfn).imag.max())

    # set the limits for the wavefunction animation
    wax.set_ylim(minwfn - 0.1 * (maxwfn - minwfn), maxwfn + 0.1 * (maxwfn - minwfn))

    # define the update function for the wavefunction animation
    update = lambda i: [
        wfnplot[j].set_ydata(wfn[i][:, j // 2].real if j % 2 == 0 else wfn[i][:, j // 2].imag) for j in range(2 * psi.shape[1])
    ]

    # make the wavefunction plot animation and set the layout
    ani = anm.FuncAnimation(wfig, update, frames=range(len(wfn)), repeat=True, interval=30); wfig.tight_layout()

# show the plot
fig.tight_layout(); plt.show()
