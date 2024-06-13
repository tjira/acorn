.PHONY: molecule

# EXAMPLE TARGETS FOR EACH MODULE ==================================================================

# target to perform an imaginary-time propagation of 1D adiabatic wavepacket
aqdyn_1d_imaginary:
	../bin/acorn_expression -d 1 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*x^2"
	../bin/acorn_expression -d 1 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-x^2)" "0"
	../bin/acorn_qdyn -d 1 -i 1000 -m 1 -o 3 -p 0 -s 0.1 --savewfn

# target to perform an imaginary-time propagation of 2D adiabatic wavepacket
aqdyn_2d_imaginary:
	../bin/acorn_expression -d 2 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x^2+y^2)"
	../bin/acorn_expression -d 2 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x^2+y^2))" "0"
	../bin/acorn_qdyn -d 2 -i 1000 -m 1 -o 3 -p 0 -s 0.1 --savewfn

# target to perform an imaginary-time propagation of 3D adiabatic wavepacket
aqdyn_3d_imaginary:
	../bin/acorn_expression -d 3 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x^2+y^2+z^2)"
	../bin/acorn_expression -d 3 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x^2+y^2+z^2))" "0"
	../bin/acorn_qdyn -d 3 -i 1000 -m 1 -o 3 -p 0 -s 0.1

# target to perform an imaginary-time propagation of 4D adiabatic wavepacket
aqdyn_4d_imaginary:
	../bin/acorn_expression -d 4 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x1^2+x2^2+x3^2+x4^2)"
	../bin/acorn_expression -d 4 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x1^2+x2^2+x3^2+x4^2))" "0"
	../bin/acorn_qdyn -d 4 -i 1000 -m 1 -o 3 -p 0 -s 0.1

# target to perform a real-time propagation of 1D adiabatic wavepacket
aqdyn_1d_real:
	../bin/acorn_expression -d 1 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*x^2"
	../bin/acorn_expression -d 1 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-x^2)" "0"
	../bin/acorn_qdyn -d 1 -i 1000 -m 1 -p 0 -s 0.1 --savewfn

# target to perform a real-time propagation of 2D adiabatic wavepacket
aqdyn_2d_real:
	../bin/acorn_expression -d 2 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x^2+y^2)"
	../bin/acorn_expression -d 2 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x^2+y^2))" "0"
	../bin/acorn_qdyn -d 2 -i 1000 -m 1 -p 0 -s 0.1 --savewfn

# target to perform an real-time propagation of 3D adiabatic wavepacket
aqdyn_3d_real:
	../bin/acorn_expression -d 3 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x^2+y^2+z^2)"
	../bin/acorn_expression -d 3 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x^2+y^2+z^2))" "0"
	../bin/acorn_qdyn -d 3 -i 1000 -m 1 -p 0 -s 0.1

# target to perform an real-time propagation of 4D adiabatic wavepacket
aqdyn_4d_real:
	../bin/acorn_expression -d 4 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x1^2+x2^2+x3^2+x4^2)"
	../bin/acorn_expression -d 4 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x1^2+x2^2+x3^2+x4^2))" "0"
	../bin/acorn_qdyn -d 4 -i 1000 -m 1 -p 0 -s 0.1

# target to perform a classical 1D real-time adiabatic propagation
cdyn_1d_adiabatic:
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

random_cdyn_1d_adiabatic_ds_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)" --adiabatic

random_cdyn_1d_diabatic_ds_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"

random_cdyn_1d_adiabatic_ts_1:
	../bin/acorn_cdyn -c -10 -e 2 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))" --adiabatic

random_cdyn_1d_diabatic_ts_1:
	../bin/acorn_cdyn -c -10 -e 1 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"

random_cdyn_1d_adiabatic_tully_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" --adiabatic

random_cdyn_1d_diabatic_tully_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 3500 -l 500 -m 2000 -p 15 -r 1 -s 1 -t 1000 -u "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"

random_qdyn_1d_ds_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_qdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic

random_qdyn_1d_ts_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "0" "0" "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_qdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic

random_qdyn_1d_tully_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_qdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic

random_dyn_1d_ds_1: random_qdyn_1d_ds_1 random_cdyn_1d_diabatic_ds_1 random_cdyn_1d_adiabatic_ds_1
random_dyn_1d_ts_1: random_qdyn_1d_ts_1 random_cdyn_1d_diabatic_ts_1 random_cdyn_1d_adiabatic_ts_1
random_dyn_1d_tully_1: random_qdyn_1d_tully_1 random_cdyn_1d_diabatic_tully_1 random_cdyn_1d_adiabatic_tully_1
