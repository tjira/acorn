#!/usr/bin/env python

import argparse as ap, itertools as it, matplotlib.animation as anm, matplotlib.pyplot as plt, numpy as np, scipy as sp, scipy.linalg

np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%20.14f" % x)); np.random.seed(0)

# EXAMPLES
# ./shcdet.py -p 10 1 -r -10 0.5 -s 1 -m 2000 -t 1 -i 5000 -v "[[0.01*np.tanh(0.6*r1),0.001*np.exp(-r1**2)],[0.001*np.exp(-r1**2),-0.01*np.tanh(0.6*r1)]]" --adiabatic -n 1000 --lzsh
# ./shcdet.py -p 10 1 -r -10 0.5 -s 2 -m 2000 -t 1 -i 5000 -v "[[0.01*np.tanh(0.6*r1),0.001*np.exp(-r1**2),0*r1],[0.001*np.exp(-r1**2),0*r1,0.001*np.exp(-r1**2)],[0*r1,0.001*np.exp(-r1**2),-0.01*np.tanh(0.6*r1)]]" --adiabatic -n 1000 --fssh
# ./shcdet.py -p 15 1 -r -15 0.5 -s 2 -m 2000 -t 1 -i 5000 -v "[[0.03*(np.tanh(1.6*r1)+np.tanh(1.6*(r1+7))),0.005*np.exp(-r1**2),0.005*np.exp(-(r1+7)**2)],[0.005*np.exp(-r1**2),-0.03*(np.tanh(1.6*r1)+np.tanh(1.6*(r1-7))),0.005*np.exp(-(r1-7)**2)],[0.005*np.exp(-(r1+7)**2),0.005*np.exp(-(r1-7)**2),-0.03*(np.tanh(1.6*(r1+7))-np.tanh(1.6*(r1-7)))]]" --adiabatic -n 1000 --fssh

# SECTION FOR PARSING COMMAND LINE ARGUMENTS ===========================================================================

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
parser.add_argument("-n", "--trajectories", help="Number of classical trajectories. (default: %(default)s)", type=int, default=1000)
parser.add_argument("-s", "--state", help="Initial state of the trajectories. (default: %(default)s)", type=int, default=0)
parser.add_argument("-t", "--timestep", help="Time step of the simulation. (default: %(default)s)", type=float, default=0.1)
parser.add_argument("-p", "--momentum", help="Momentum distribution. (default: %(default)s)", type=float, nargs= "+", default=[0, 1])
parser.add_argument("-r", "--position", help="Coords distribution. (default: %(default)s)", type=float, nargs= "+", default=[1, 0.5])
parser.add_argument("-u", "--adiabatic", help="Adiabatic transform. (default: %(default)s)", action=ap.BooleanOptionalAction)
parser.add_argument("-v", "--potential", help="Potential matrix. (default: %(default)s)", type=str, default="[[0.5*r1**2]]")

# surface hopping arguments
parser.add_argument("--lzsh", help="Enable LZSH algorithm. (default: %(default)s)", action=ap.BooleanOptionalAction)
parser.add_argument("--fssh", help="Enable FSSH algorithm. (default: %(default)s)", action=ap.BooleanOptionalAction)
parser.add_argument("--mash", help="Enable MASH algorithm. (default: %(default)s)", action=ap.BooleanOptionalAction)

# parse the arguments
args = parser.parse_args()

# replace the variables in the string expressions
for i in range(1, 10): args.potential = args.potential.replace(f"r{i}", f"r[:, {i - 1}]")

# print the help message if the flag is set
if args.help: parser.print_help(); exit()

# SURFACE HOPPING ALGORITHMS =====================================================================================================

# define the state vector and containers for the potential and transformation matrices
s = np.zeros((args.trajectories, args.iterations + 1), dtype=int) + args.state; Vs, Us = [], []

# define the Fewest Switches Surface Hopping algorithm
def fssh(i, substeps=10):

    # fill the coefficients on the first iteration
    if i == 1: 
        for j in range(args.trajectories): C[j, s[j, 0]] = 1

    # compute the eigenvector overlap and the trajectory indices
    O = Us[-1].swapaxes(1, 2) @ Us[-2]; traji = np.arange(args.trajectories)

    # calculate the time derivative coupling
    TDC = (np.transpose(O, (0, 2, 1)) - O) / (2 * args.timestep)

    # vectorized derivative function
    dC = lambda C: -1j * np.einsum("ijj->ij", Vs[-1]) * C - np.einsum("ijk,ik->ij", TDC, C)

    # loop over the substeps
    for _ in range(substeps):

        # calculate the Runge-Kutta coefficients
        k1 = dC(C                                      )
        k2 = dC(C + 0.5 * args.timestep * k1 / substeps)
        k3 = dC(C + 0.5 * args.timestep * k2 / substeps)
        k4 = dC(C + 1.0 * args.timestep * k3 / substeps)

        # update the electronic coefficients and get trajectory indices
        C[:] = C + args.timestep * (k1 + 2 * k2 + 2 * k3 + k4) / substeps / 6.0;

        # calculate the denominator for the hopping probability
        denominator = 0.5 * (np.abs(C[traji, s[:, i]])**2 + 1e-14) * substeps / args.timestep

        # calculate the hopping probability
        p = TDC[traji, s[:, i], :] * (C * np.conjugate(C[traji, s[:, i]][:, None])).real / denominator[:, None]

        # calculate the hopping mask
        hopmask = (np.random.rand(args.trajectories)[:, None] < p)

        # perform the hop using the mask to the state that has maximum hopping probability
        if hopmask.any(): s[hopmask.any(axis=1), i:] = np.argmax(p, axis=1)[hopmask.any(axis=1), None]

# define the Landau-Zener Surface Hopping algorithm
def lzsh(i):

    # loop over the trajectories
    for j in range(args.trajectories):

        # generate random number
        rn = np.random.rand()

        # loop over the states
        for k in (l for l in range(V.shape[1]) if l != s[j, i]):

            # skip current state
            if k == s[j, i]: continue

            # calculate the energy differences
            Z0 = abs(Vs[-1][j, k, k] - Vs[-1][j, s[j, i], s[j, i]])
            Z1 = abs(Vs[-2][j, k, k] - Vs[-2][j, s[j, i], s[j, i]])
            Z2 = abs(Vs[-3][j, k, k] - Vs[-3][j, s[j, i], s[j, i]])

            # calculate first derivatives of the energy differences
            dZ0, dZ1 = (Z0 - Z1) / args.timestep, (Z1 - Z2) / args.timestep

            # calculate the second derivative
            ddZ0 = (Z0 - 2 * Z1 + Z2) / args.timestep**2

            # check if the trajectory is in the place for a hop
            if (dZ0 * dZ1 > 0 or (dZ0 * dZ1 < 0 and ddZ0 < 0)): continue

            # calculate the hopping probability
            p = np.exp(-0.5 * np.pi * np.sqrt(Z0**3 / ddZ0)) if (Z0**3 * ddZ0 > 0) else 0

            # hop to another state
            if rn < p:
                s[j, i:] = k; break

# define the Mapping Approach to Surface Hopping algorithm
def mash(i):

    # sample the bloch vectors if i == 1
    for j in (k for k in range(args.trajectories) if i == 1):

        # define the angles
        p, ct = 2 * np.pi * np.random.rand(), np.sqrt(np.random.rand()); st = np.sqrt(1 - ct**2)

        # fill the row
        for l in range(len(SPHR)): SPHR[l, j] = np.array([st * np.cos(p), st * np.sin(p), ct])

    # compute the eigenvector overlap and store spheres
    O = Us[-1].swapaxes(1, 2) @ Us[-2]; SPHRP = SPHR.copy()

    # calculate the time derivative coupling
    TDC = (np.transpose(O, (0, 2, 1)) - O) / (2 * args.timestep)

    # loop over the trajectories
    for j in range(args.trajectories):

        # loop over all spheres
        for k in range(len(SPHR)):

            # get upper and lower index
            su, sd = max(SIND[j, k]), min(SIND[j, k])

            # calculate the propagator
            omega = sp.linalg.expm(np.array([
                [0,                                     Vs[-1][j, sd, sd] - Vs[-1][j, su, su], 2 * TDC[j, sd, su]],
                [Vs[-1][j, su, su] - Vs[-1][j, sd, sd], 0,                                     0                 ],
                [-2 * TDC[j, sd, su],                   0,                                     0                 ]
            ]) * args.timestep)

            # propagate the Bloch vector
            SPHR[k, j] = omega @ SPHR[k, j]

        # loop over spheres
        for k in range(len(SPHR)):

            # perform the jump if the sphere crossed the equator
            if SPHRP[k, j, 2] * SPHR[k, j, 2] < 0: s[j, i:] = SIND[j, k][0] if SIND[j, k][1] == s[j, i] else SIND[j, k][1]

        # loop over sphere indices if the jump occured
        for k in (l for l in range(len(SIND[j])) if s[j, i] != s[j, i - 1]):

            # update the indices of the spheres
            if SIND[j, k][0] == s[j, i - 1] and SIND[j, k][1] != s[j, i]: SIND[j, k][0] = s[j, i]
            if SIND[j, k][1] == s[j, i - 1] and SIND[j, k][0] != s[j, i]: SIND[j, k][1] = s[j, i]

# PERFORM THE CLASSICAL DYNAMICS =================================================================================================

# create the initial conditions for the trajectories
r = np.column_stack([np.random.normal(mu, sigma, len(s)) for mu, sigma in zip(args.position[0::2], args.position[1::2])])
v = np.column_stack([np.random.normal(mu, sigma, len(s)) for mu, sigma in zip(args.momentum[0::2], args.momentum[1::2])])

# extract the number of states
nstate = np.array(eval(args.potential)).shape[0]

# create the electronic coefficients and one Bloch sphere
C = np.zeros((args.trajectories, nstate), dtype=complex)
S = np.zeros((args.trajectories, 3     ), dtype=float  )

# create the array of spheres and sphere indices for the MASH algorithm
SPHR = np.array([ S.copy()     for i in range(nstate - 1)                ]                                   )
SIND = np.array([[[s[0, 0], i] for i in range(nstate - 0) if i != s[0, 0]] for j in range(args.trajectories)])

# divide momentum by mass and define acceleration
v /= args.mass; a = np.zeros_like(r); diff = 0.0001;

# print the propagation header
print("%6s %12s %12s %12s" % ("ITER", "EKIN", "EPOT", "ETOT"))

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

        Vplusi = np.einsum("ijj->ij", Vplus)[np.arange(args.trajectories), s[:, i]]
        Vminsi = np.einsum("ijj->ij", Vmins)[np.arange(args.trajectories), s[:, i]]

        # calculate the acceleration
        a[:, j] = 0.5 * (Vminsi - Vplusi) / (diff * args.mass)

        # update the positions and velocities
        v[:, j] += 0.5 * (a[:, j] + ap[:, j]) * args.timestep; r[:, j] += (v[:, j] + 0.5 * a[:, j] * args.timestep) * args.timestep

    # calculate the diabatic potential
    V = np.array(eval(args.potential)).transpose(2, 0, 1)

    # adiabatize the potential
    E, U = np.linalg.eigh(V); V = np.einsum("ij,jk->ijk", E, np.eye(V.shape[1])); Vs.append(V); Us.append(U)

    # maximize overlap between the current and previous eigenvectors
    if i > 1: Us[-1] = Us[-1] * np.where(np.sum(Us[-1] * Us[-2], axis=1) < 0, -1, 1)[:, None, :]

    # surface hopping algorithm
    if (args.fssh and i > 0): fssh(i)
    if (args.lzsh and i > 1): lzsh(i)
    if (args.mash and i > 0): mash(i)

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

# PRINT AND PLOT THE RESULTS =====================================================================================================

# print the final populations
print(); print("FINAL POPULATION: %s" % (np.array([(s[:, -1] == j).sum() / s.shape[0] for j in range(V.shape[1])])))

# create the subplots
fig, axs = plt.subplots(1, 1, figsize=(8, 6));

# plot the population
axs.plot(np.arange(args.iterations + 1) * args.timestep, [[(s[:, i] == j).sum() / s.shape[0] for j in range(V.shape[1])] for i in range(s.shape[1])])

# set the labels
axs.set_xlabel("Time (a.u.)"); axs.set_ylabel("Population")

# show the plot
plt.tight_layout(); plt.show()
