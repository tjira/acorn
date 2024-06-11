.PHONY: molecule

# EXAMPLE TARGETS FOR EACH MODULE ==================================================================

# target to generate the 1D adiabatic wavepacket guess
aqdyn_guess:
	../bin/acorn_expression -g -16 16 -o PSI_ADIA_GUESS.mat -p 1024 -e "exp(-(x-1)^2)" "0"

# target to perform an imaginary-time propagation of 1D adiabatic wavepacket
aqdyn_imaginary: aqdyn_guess aqdyn_potential
	../bin/acorn_aqdyn -i 200 -m 1 -n 3 -p 0 -s 0.1 --imaginary --savewfn

# target to generate the 1D adiabatic potential matrix
aqdyn_potential:
	../bin/acorn_expression -g -16 16 -o U_ADIA.mat -p 1024 -e "0.5*x^2"

# target to perform a real-time propagation of 1D adiabatic wavepacket
aqdyn_real: aqdyn_guess aqdyn_potential
	../bin/acorn_aqdyn -i 1000 -m 1 -n 1 -p 0 -s 0.1 --savewfn

# target to perform a classical 1D real-time adiabatic propagation
cdyn_adiabatic:
	../bin/acorn_cdyn -c 3 -e 0 -i 1000 -m 1 -p 0 -r 1 -s 0.01 -t 5 -u "0.5*x^2" --savetraj

# target to perform a MP3 calculation
fci: transform
	../bin/acorn_ci

# target to perform a Hartree-Fock calculation
hf: integral
	../bin/acorn_hf -d 5 -i 100 -t 1e-8

# target to calculate the integrals in specified basis
integral: molecule
	../bin/acorn_integral -b sto-3g

# target to copy the molecule file to the current folder
molecule:
	cp molecule/water.xyz molecule.xyz

# target to perform a MP2 calculation
mp2: transform
	../bin/acorn_mp -o 2

# target to perform a MP3 calculation
mp3: transform
	../bin/acorn_mp -o 3

# target to transform atomic integrals to the chosen basis
transform: hf
	../bin/acorn_transform -s

# RANDOM TARGETS ===================================================================================

random_cdyn_adiabatic_ds_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)" --adiabatic

random_cdyn_diabatic_ds_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"

random_cdyn_adiabatic_ts_1:
	../bin/acorn_cdyn -c -10 -e 2 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))" --adiabatic

random_cdyn_diabatic_ts_1:
	../bin/acorn_cdyn -c -10 -e 1 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"

random_cdyn_adiabatic_tully_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" --adiabatic

random_cdyn_diabatic_tully_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"

random_nqdyn_ds_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_nqdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic

random_nqdyn_ts_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "0" "0" "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_nqdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic

random_nqdyn_tully_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_nqdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic

random_dyn_ds_1: random_nqdyn_ds_1 random_cdyn_diabatic_ds_1 random_cdyn_adiabatic_ds_1
random_dyn_ts_1: random_nqdyn_ts_1 random_cdyn_diabatic_ts_1 random_cdyn_adiabatic_ts_1
random_dyn_tully_1: random_nqdyn_tully_1 random_cdyn_diabatic_tully_1 random_cdyn_adiabatic_tully_1
