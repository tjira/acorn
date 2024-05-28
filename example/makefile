.PHONY: molecule integral

# target to generate the 1D adiabatic wavepacket guess
aqdyn_guess:
	acorn_expression -g -16 16 -o PSI_ADIA_GUESS.mat -p 2048 -e "exp(-(x-1)^2)" "0"

# target to perform an imaginary-time propagation of 1D adiabatic wavepacket
aqdyn_imaginary: aqdyn_guess aqdyn_potential
	acorn_aqdyn -i 200 -m 1 -n 3 -p 0 -s 0.1 --imaginary --savewfn

# target to generate the 1D adiabatic potential matrix
aqdyn_potential:
	acorn_expression -g -16 16 -o U_ADIA.mat -p 2048 -e "0.5*x^2"

# target to perform a real-time propagation of 1D adiabatic wavepacket
aqdyn_real: aqdyn_guess aqdyn_potential
	acorn_aqdyn -i 1000 -m 1 -n 1 -p 0 -s 0.1 --savewfn

# target to perform a classical 1D real-time adiabatic propagation
cdyn_adiabatic:
	acorn_cdyn -c 3 -i 1000 -m 1 -p 0 -r 1 -s 0.01 -t 5 -u "0.5*x^2" -d "x" --savetraj

# target to perform a classical 1D real-time nonadiabatic propagation
cdyn_nonadiabatic:
	acorn_cdyn -c -10 -i 350 -m 2000 -p 10.95 -r 1 -s 10 -t 100 -u "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)" -d "0.006/cosh(0.6*x)^2" "(-0.002)*x*exp(-x^2)" "(-0.002)*x*exp(-x^2)" "(-0.006)/cosh(0.6*x)^2" --savetraj

# target to perform a Hartree-Fock calculation
hf: integral
	acorn_hf -i 100 -t 1e-8

# target to calculate the integrals in specified basis
integral: molecule
	acorn_integral -b sto-3g

# target to copy the molecule file to the current folder
molecule:
	cp molecule/water.xyz molecule.xyz

# target to perform a MP2 calculation
mp2: hf
	acorn_mp -o 2

# target to perform a MP3 calculation
mp3: hf
	acorn_mp -o 3

# target to generate the 1D diabatic wavepacket guess
nqdyn_guess:
	acorn_expression -g -24 24 -o PSI_DIA_GUESS.mat -p 4096 -e "exp(-(x+10)^2)" "0" "0" "0"

# target to generate the 1D diabatic potential matrix
nqdyn_potential:
	acorn_expression -g -24 24 -o U_DIA.mat -p 4096 -e "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"

# target to perform a real-time propagation of 1D diabatic wavepacket
nqdyn_real: nqdyn_guess nqdyn_potential
	acorn_nqdyn -i 350 -m 2000 -p 10.95 -s 10 --adiabatic --savewfn

# target to transform atomic integrals to the chosen basis
transform: hf
	acorn_transform -s
