# import plotting and animating functions
import matplotlib.animation as anm
import matplotlib.pyplot as plt

# import numpy library
import numpy as np

ngrid, dt, iters, m = 512, 32, 512, 2048

def V(x):
    # define the matrix
    V = np.zeros([2, 2])

    # define the potential
    V[0,0] = 0.01 * np.tanh(0.6 * x); V[1,1] = -V[0,0]
    V[0,1] = 0.001 * np.exp(-x**2); V[1,0] = V[0,1]

    # return V
    return V

if __name__ == "__main__":
    # define the grid and the grid spacing
    r = np.linspace(-32, 32, ngrid); dr = r[1] - r[0]

    # define the kinetic energy matrix, Hamiltonian, potential in adiabatic basis and the transformation matrix
    T, H, VV, UU = np.zeros([ngrid, ngrid]), np.zeros([2 * ngrid, 2 * ngrid]), list(), np.zeros((2, 2, ngrid))

    # fill the kinetic energy matrix
    for i in range(ngrid):
        for j in range(ngrid):
            T[i,j] = np.pi ** 2 / 3 if i == j else 2 / (i - j) ** 2
            T[i,j] *= 1 / (2 * m * dr ** 2) * (-1) ** (i - j)

    # assign the kinetic energy matrix to the Hamiltonian
    H[:ngrid,:ngrid] = T; H[ngrid:2 * ngrid, ngrid:2 * ngrid] = T

    for i in range(ngrid):
        # add the potential energy
        H[i, i + ngrid] += V(r[i])[0, 1]; H[i + ngrid, i] += V(r[i])[1, 0]
        H[i, i] += V(r[i])[0, 0]; H[i + ngrid, i + ngrid] += V(r[i])[1, 1]

        # find the adiabatic basis
        eigv, eigf = np.linalg.eigh(V(r[i]))

        # assign the adiabatic results
        UU[:, :, i] = eigf; VV.append(eigv)

    # solve the time independent Schrodinger equation
    eigv, eigf = np.linalg.eigh(H)

    # define the initial wavefunction and an array of adiabatic wfns
    psi, psis = np.zeros(2 * ngrid, dtype=complex), []

    # fill the initial wavefunction in diabatic basis
    psi[:ngrid] = np.exp(-(r + 15)**2) * np.exp(1j * np.sqrt(0.06 * m) * r); psi /= np.sqrt(np.vdot(psi, psi) * dr)

    # define the array for density matrices
    rho = np.zeros((2, 2, iters), dtype=complex)

    # calculate the expansion coefficients of the diabatic wavefunction
    ci = np.array([sum(psi * eigf[:, i]) * dr for i in range(2 * ngrid)])

    max_y = np.max(abs(psi))

    for i in range(iters):
        # propagate the wavefunction in the diabatic basis
        psi = np.array(sum([ci[j] * eigf[:, j] * np.exp(-1j * eigv[j] * i * dt) for j in range(2 * ngrid)]))

        # normalize the wavefunction
        psi /= np.sqrt(np.vdot(psi, psi) * dr)

        # extract the adiabatic basis
        chi_ad = np.array([np.transpose(UU[:,:,i])@[psi[i], psi[i + ngrid]] for i in range(ngrid)])

        # calculate the density matrix
        rho[:, :, i] = sum([np.array([[abs(chi_ad[i, 0])**2, np.conj(chi_ad[i, 1]) * chi_ad[i, 0]],[np.conj(chi_ad[i, 0]) * chi_ad[i, 1], abs(chi_ad[i, 1])**2]]) for i in range(ngrid)])

        # append the wavefunction in the adiabatic basis
        psis.append(chi_ad.T.reshape(2 * ngrid))

    # define the animation plot parameters
    [fig, ax] = plt.subplots(2, 2); plt.tight_layout()

    # plot the adiabatic potentials
    ax[0, 1].plot(r, np.array(VV)[:, 0])
    ax[0, 1].plot(r, np.array(VV)[:, 1])

    # plot the density matrix elements
    ax[1, 1].plot(np.arange(iters) * dt, abs(rho[0, 0, :]), label=r'$\rho_{el}^{00}$')
    ax[1, 1].plot(np.arange(iters) * dt, abs(rho[1, 1, :]), label=r'$\rho_{el}^{11}$')
    ax[1, 1].plot(np.arange(iters) * dt, abs(rho[0, 1, :]), label=r'$|\rho_{el}^{01}|$')
    ax[1, 1].legend()

    # set limits of all y axis
    ax[0, 0].set_ylim(0, max([np.max(abs(psi)) for psi in psis]))
    ax[1, 0].set_ylim(0, max([np.max(abs(psi)) for psi in psis]))

    # define the plot objects
    psiplots = [ax[0, 0].plot(r, abs(psis[0][ngrid:]))[0], ax[1, 0].plot(r, abs(psis[0][:ngrid]))[0]]

    # define the update function
    def update(j):
        psiplots[0].set_ydata(abs(psis[j][ngrid:]))
        psiplots[1].set_ydata(abs(psis[j][:ngrid]))

    # animate the wavefunction
    ani = anm.FuncAnimation(fig, update, frames=len(psis), interval=30); plt.show(); plt.close("all")
