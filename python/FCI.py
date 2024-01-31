"""
A Psi4 input script to compute Full Configuration Interaction from a SCF reference

Requirements:
SciPy 0.13.0+, NumPy 1.7.2+

References:
Equations from [Szabo:1996]
"""

__authors__ = "Tianyuan Zhang"
__credits__ = ["Tianyuan Zhang", "Jeffrey B. Schriber", "Daniel G. A. Smith"]

__copyright__ = "(c) 2014-2018, The Psi4NumPy Developers"
__license__ = "BSD-3-Clause"
__date__ = "2017-05-26"

import time
import numpy as np
np.set_printoptions(precision=5, linewidth=200, suppress=True)
import psi4

# Check energy against psi4?
compare_psi4 = True

# Memory for Psi4 in GB
# psi4.core.set_memory(int(2e9), False)
psi4.core.set_output_file('output.dat', False)

# Memory for numpy in GB
numpy_memory = 2

mol = psi4.geometry("""
 H  0.14388500725555  0.00000000000000  0.00000000000000
 H  0.85611499274445  0.00000000000000  0.00000000000000
symmetry c1
""")

# mol = psi4.geometry("""
#  O -0.04258229224928 -0.01476311599192 -0.01697740475215
#  H  0.65343210274828 -0.51640346349129  0.47582054313492
#  H  0.38535018950100  0.87626657948321 -0.06024313838278
# symmetry c1
# """)

mol = psi4.geometry("""
 H  0.02226951616055  0.00000000000000  0.00000000000000
 F  0.97773048383945  0.00000000000000  0.00000000000000
symmetry c1
""")


psi4.set_options({'basis': 'sto-3g',
                  'scf_type': 'pk',
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8})

print('\nStarting SCF and integral build...')
t = time.time()

# First compute SCF energy using Psi4
scf_e, wfn = psi4.energy('SCF', return_wfn=True)

# Grab data from wavfunction class
# C = wfn.Ca()
C = psi4.core.Matrix.from_array(np.loadtxt("example/input/C.mat"))
ndocc = wfn.doccpi()[0]
nmo = wfn.nmo()

# Compute size of Hamiltonian in GB
from scipy.special import comb
nDet = comb(nmo, ndocc)**2
H_Size = nDet**2 * 8e-9
print('\nSize of the Hamiltonian Matrix will be %4.2f GB.' % H_Size)
if H_Size > numpy_memory:
    clean()
    raise Exception("Estimated memory utilization (%4.2f GB) exceeds numpy_memory \
                    limit of %4.2f GB." % (H_Size, numpy_memory))

# Integral generation from Psi4's MintsHelper
t = time.time()
mints = psi4.core.MintsHelper(wfn.basisset())
# H = np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())
H = np.loadtxt("example/input/T.mat") + np.loadtxt("example/input/V.mat")

print('\nTotal time taken for ERI integrals: %.3f seconds.\n' % (time.time() - t))

#Make spin-orbital MO
print('Starting AO -> spin-orbital MO transformation...')
t = time.time()
MO = np.asarray(mints.mo_spin_eri(C, C))
MOO = np.asarray(mints.mo_eri(C, C, C, C))

# Update H, transform to MO basis and tile for alpha/beta spin
H = np.einsum('uj,vi,uv', C, C, H)
H = np.repeat(H, 2, axis=0)
H = np.repeat(H, 2, axis=1)

# Make H block diagonal
spin_ind = np.arange(H.shape[0], dtype=int) % 2
H *= (spin_ind.reshape(-1, 1) == spin_ind)

print("Hms", ((H - np.loadtxt("example/input/HMS.mat"))**2).sum())
print("Jms", ((MO.reshape(MO.shape[0]**2, MO.shape[1]**2) - np.loadtxt("example/input/JMS.mat"))**2).sum())

# print(C)
# print(np.loadtxt("JMS.mat"))
# print(MO.reshape(16, 16))

print('..finished transformation in %.3f seconds.\n' % (time.time() - t))

from helper_CI import Determinant, HamiltonianGenerator
from itertools import combinations

print('Generating %d Full CI Determinants...' % (nDet))
t = time.time()
detList = []
for alpha in combinations(range(nmo), ndocc):
    for beta in combinations(range(nmo), ndocc):
        detList.append(Determinant(alphaObtList=alpha, betaObtList=beta))

print('..finished generating determinants in %.3f seconds.\n' % (time.time() - t))

print('Generating Hamiltonian Matrix...')

t = time.time()
Hamiltonian_generator = HamiltonianGenerator(H, MO)
Hamiltonian_matrix = Hamiltonian_generator.generateMatrix(detList)

print(Hamiltonian_matrix)
print(np.loadtxt("example/input/HCI.mat"))
print((Hamiltonian_matrix - np.loadtxt("example/input/HCI.mat")).sum())

print('..finished generating Matrix in %.3f seconds.\n' % (time.time() - t))

print('Diagonalizing Hamiltonian Matrix...')

t = time.time()

e_fci, wavefunctions = np.linalg.eigh(Hamiltonian_matrix)
print('..finished diagonalization in %.3f seconds.\n' % (time.time() - t))

fci_mol_e = e_fci[0] + mol.nuclear_repulsion_energy()

print('# Determinants:     % 16d' % (len(detList)))

print('SCF energy:         % 16.10f' % (scf_e))
print('FCI correlation:    % 16.10f' % (fci_mol_e - scf_e))
print('Total FCI energy:   % 16.10f' % (fci_mol_e))

if compare_psi4:
    psi4.compare_values(psi4.energy('FCI'), fci_mol_e, 6, 'FCI Energy')
