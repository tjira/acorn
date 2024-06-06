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
	../bin/acorn_cdyn -c -10 -e 0 -i 350 -m 2000 -p 15 -r 1 -s 10 -t 1000 -u "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)" --adiabatic
	gnuplot -e "columns=2; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_LZ_ADIA'" -c ../script/plot.gp P_LZ_ADIA.mat

random_cdyn_diabatic_ds_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 350 -m 2000 -p 15 -r 1 -s 10 -t 1000 -u "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"
	gnuplot -e "columns=2; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_LZ_DIA'" -c ../script/plot.gp P_LZ_DIA.mat

random_cdyn_adiabatic_ts_1:
	../bin/acorn_cdyn -c -10 -e 2 -i 350 -m 2000 -p 15 -r 1 -s 10 -t 1000 -u "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))" --adiabatic
	gnuplot -e "columns=3; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_LZ_ADIA'" -c ../script/plot.gp P_LZ_ADIA.mat

random_cdyn_diabatic_ts_1:
	../bin/acorn_cdyn -c -10 -e 1 -i 350 -m 2000 -p 15 -r 1 -s 10 -t 1000 -u "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"
	gnuplot -e "columns=3; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_LZ_DIA'" -c ../script/plot.gp P_LZ_DIA.mat

random_cdyn_adiabatic_tully_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 350 -m 2000 -p 15 -r 1 -s 10 -t 1000 -u "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" --adiabatic
	gnuplot -e "columns=2; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_LZ_ADIA'" -c ../script/plot.gp P_LZ_ADIA.mat

random_cdyn_diabatic_tully_1:
	../bin/acorn_cdyn -c -10 -e 0 -i 350 -m 2000 -p 15 -r 1 -s 10 -t 1000 -u "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"
	gnuplot -e "columns=2; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_LZ_DIA'" -c ../script/plot.gp P_LZ_DIA.mat

random_nqdyn_ds_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_nqdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic
	gnuplot -e "columns=4; only='1 4'; xmin=-20; xmax=20; ymin=-0.011; ymax=0.011; format='png'; output='U_DIA'" -c ../script/plot.gp U_DIA.mat
	gnuplot -e "columns=4; only='1 4'; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_DIA'" -c ../script/plot.gp P_DIA.mat
	gnuplot -e "columns=2; xmin=-20; xmax=20; ymin=-0.011; ymax=0.011; format='png'; output='U_ADIA'" -c ../script/plot.gp U_ADIA.mat
	gnuplot -e "columns=4; only='1 4'; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_ADIA'" -c ../script/plot.gp P_ADIA.mat

random_nqdyn_ts_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-x^2)" "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)" "0.005*exp(-(x+7)^2)" "0.005*exp(-(x-7)^2)" "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "0" "0" "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_nqdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic
	gnuplot -e "columns=9; only='1 5 9'; xmin=-20; xmax=20; ymin=-0.07; ymax=0.07; format='png'; output='U_DIA'" -c ../script/plot.gp U_DIA.mat
	gnuplot -e "columns=9; only='1 5 9'; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_DIA'" -c ../script/plot.gp P_DIA.mat
	gnuplot -e "columns=3; xmin=-20; xmax=20; ymin=-0.07; ymax=0.07; format='png'; output='U_ADIA'" -c ../script/plot.gp U_ADIA.mat
	gnuplot -e "columns=9; only='1 5 9'; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_ADIA'" -c ../script/plot.gp P_ADIA.mat

random_nqdyn_tully_1:
	../bin/acorn_expression -g -24 40 -o U_DIA.mat -p 1024 -e "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)" "0.005*exp(-x^2)" "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"
	../bin/acorn_expression -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	../bin/acorn_nqdyn -i 350 -m 2000 -p 15 -s 10 --adiabatic
	gnuplot -e "columns=4; only='1 4'; xmin=-20; xmax=20; ymin=-0.011; ymax=0.011; format='png'; output='U_DIA'" -c ../script/plot.gp U_DIA.mat
	gnuplot -e "columns=4; only='1 4'; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_DIA'" -c ../script/plot.gp P_DIA.mat
	gnuplot -e "columns=2; xmin=-20; xmax=20; ymin=-0.011; ymax=0.011; format='png'; output='U_ADIA'" -c ../script/plot.gp U_ADIA.mat
	gnuplot -e "columns=4; only='1 4'; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; format='png'; output='P_ADIA'" -c ../script/plot.gp P_ADIA.mat

random_dyn_ds_1: random_nqdyn_ds_1 random_cdyn_diabatic_ds_1 random_cdyn_adiabatic_ds_1
	awk '{print $$1, $$2, $$5}' P_DIA.mat > P_DIA_PO.mat && awk '{print $$1, $$2, $$5}' P_ADIA.mat > P_ADIA_PO.mat
	gnuplot -e "columns=2; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; titles='EX1 EX2 EX3 LZ1 LZ2 LZ3'; format='png'; output='P_MIX_ADIA'" -c ../script/plot.gp P_ADIA_PO.mat P_LZ_ADIA.mat
	gnuplot -e "columns=2; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; titles='EX1 EX2 EX3 LZ1 LZ2 LZ3'; format='png'; output='P_MIX_DIA'" -c ../script/plot.gp P_DIA_PO.mat P_LZ_DIA.mat

random_dyn_ts_1: random_nqdyn_ts_1 random_cdyn_diabatic_ts_1 random_cdyn_adiabatic_ts_1
	awk '{print $$1, $$2, $$6, $$10}' P_DIA.mat > P_DIA_PO.mat && awk '{print $$1, $$2, $$6, $$10}' P_ADIA.mat > P_ADIA_PO.mat
	gnuplot -e "columns=3; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; titles='EX1 EX2 EX3 LZ1 LZ2 LZ3'; format='png'; output='P_MIX_ADIA'" -c ../script/plot.gp P_ADIA_PO.mat P_LZ_ADIA.mat
	gnuplot -e "columns=3; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; titles='EX1 EX2 EX3 LZ1 LZ2 LZ3'; format='png'; output='P_MIX_DIA'" -c ../script/plot.gp P_DIA_PO.mat P_LZ_DIA.mat

random_dyn_tully_1: random_nqdyn_tully_1 random_cdyn_diabatic_tully_1 random_cdyn_adiabatic_tully_1
	awk '{print $$1, $$2, $$5}' P_DIA.mat > P_DIA_PO.mat && awk '{print $$1, $$2, $$5}' P_ADIA.mat > P_ADIA_PO.mat
	gnuplot -e "columns=2; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; titles='EX1 EX2 LZ1 LZ2'; format='png'; output='P_MIX_ADIA'" -c ../script/plot.gp P_ADIA_PO.mat P_LZ_ADIA.mat
	gnuplot -e "columns=2; xmin=0; xmax=3500; ymin=-0.05; ymax=1.05; titles='EX1 EX2 LZ1 LZ2'; format='png'; output='P_MIX_DIA'" -c ../script/plot.gp P_DIA_PO.mat P_LZ_DIA.mat
