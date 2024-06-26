.PHONY: molecule

# target to clean the folder
clean:
	rm -f *.mat *.xyz

# target to copy the molecule file to the current folder
molecule:
	cp "molecule/water.xyz" molecule.xyz

# EXAMPLE TARGETS FOR EACH MODULE ==================================================================

# target to perform a classical 1D real-time propagation
cdyn_1d_HO:
	"../bin/acorn_cdyn" -c 3 -e 0 -i 1000 -m 1 -p 0 -r 1 -s 0.01 -t 5 -u "0.5*x^2" --savetraj

# target to perform an imaginary-time propagation of 1D adiabatic wavepacket
qdyn_1d_HO_imaginary:
	"../bin/acorn_expression" -d 1 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*x^2"
	"../bin/acorn_expression" -d 1 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-x^2)" "0"
	"../bin/acorn_qdyn" -d 1 -i 1000 -m 1 -o 3 -p 0 -s 0.1 --savewfn

# target to perform an imaginary-time propagation of 2D adiabatic wavepacket
qdyn_2d_HO_imaginary:
	"../bin/acorn_expression" -d 2 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x^2+y^2)"
	"../bin/acorn_expression" -d 2 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x^2+y^2))" "0"
	"../bin/acorn_qdyn" -d 2 -i 1000 -m 1 -o 3 -p 0 -s 0.1 --savewfn

# target to perform an imaginary-time propagation of 3D adiabatic wavepacket
qdyn_3d_HO_imaginary:
	"../bin/acorn_expression" -d 3 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x^2+y^2+z^2)"
	"../bin/acorn_expression" -d 3 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x^2+y^2+z^2))" "0"
	"../bin/acorn_qdyn" -d 3 -i 1000 -m 1 -o 3 -p 0 -s 0.1

# target to perform an imaginary-time propagation of 4D adiabatic wavepacket
qdyn_4d_HO_imaginary:
	"../bin/acorn_expression" -d 4 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x1^2+x2^2+x3^2+x4^2)"
	"../bin/acorn_expression" -d 4 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x1^2+x2^2+x3^2+x4^2))" "0"
	"../bin/acorn_qdyn" -d 4 -i 1000 -m 1 -o 3 -p 0 -s 0.1

# target to perform a real-time propagation of 1D adiabatic wavepacket
qdyn_1d_HO_real:
	"../bin/acorn_expression" -d 1 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*x^2"
	"../bin/acorn_expression" -d 1 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-(x-1)^2)" "0"
	"../bin/acorn_qdyn" -d 1 -i 1000 -m 1 -p 0 -s 0.1 --savewfn

# target to perform a real-time propagation of 2D adiabatic wavepacket
qdyn_2d_HO_real:
	"../bin/acorn_expression" -d 2 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x^2+y^2)"
	"../bin/acorn_expression" -d 2 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-((x-1)^2+(y-1)^2))" "0"
	"../bin/acorn_qdyn" -d 2 -i 1000 -m 1 -p 0 -s 0.1 --savewfn

# target to perform an real-time propagation of 3D adiabatic wavepacket
qdyn_3d_HO_real:
	"../bin/acorn_expression" -d 3 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x^2+y^2+z^2)"
	"../bin/acorn_expression" -d 3 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-((x-1)^2+(y-1)^2+(z-1)^2))" "0"
	"../bin/acorn_qdyn" -d 3 -i 1000 -m 1 -p 0 -s 0.1

# target to perform an real-time propagation of 4D adiabatic wavepacket
qdyn_4d_HO_real:
	"../bin/acorn_expression" -d 4 -g -8 8 -o U_DIA.mat -p 64 -e "0.5*(x1^2+x2^2+x3^2+x4^2)"
	"../bin/acorn_expression" -d 4 -g -8 8 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-((x1-1)^2+(x2-1)^2+(x3-1)^2+(x4-1)^2))" "0"
	"../bin/acorn_qdyn" -d 4 -i 1000 -m 1 -p 0 -s 0.1

# target to perform a Hartree-Fock calculation
hf: molecule
	"../bin/acorn_integral" -b sto-3g && ../bin/acorn_hf -d 5 -i 100 -t 1e-8

# target to perform a MP2 calculation
mp2: molecule
	"../bin/acorn_integral" -b sto-3g && ../bin/acorn_hf -d 5 -i 100 -t 1e-8 && ../bin/acorn_transform -s && ../bin/acorn_mp -o 2

# target to perform a MP3 calculation
mp3: molecule
	"../bin/acorn_integral" -b sto-3g && ../bin/acorn_hf -d 5 -i 100 -t 1e-8 && ../bin/acorn_transform -s && ../bin/acorn_mp -o 3

# target to perform a FCI calculation
fci: molecule
	"../bin/acorn_integral" -b sto-3g && ../bin/acorn_hf -d 5 -i 100 -t 1e-8 && ../bin/acorn_transform -s && ../bin/acorn_ci

# RANDOM TARGETS ===================================================================================

TULLY_1 = "if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))" "0.005*exp(-x^2)"                                        \
          "0.005*exp(-x^2)"                                   "(-1)*if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"

TULLY_2 = "0"                    "0.015*exp(-0.06*x^2)"       \
          "0.015*exp(-0.06*x^2)" "(-0.1)*exp(-0.28*x^2)+0.05"

DS_1 = "0.01*tanh(0.6*x)" "0.001*exp(-x^2)"     \
       "0.001*exp(-x^2)"  "(-0.01)*tanh(0.6*x)"

TS_1 = "0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)"                       "0.005*exp(-(x+7)^2)"                       \
       "0.005*exp(-x^2)"                    "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7)))" "0.005*exp(-(x-7)^2)"                       \
       "0.005*exp(-(x+7)^2)"                "0.005*exp(-(x-7)^2)"                   "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"

TS_2 = "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.005*exp(-x^2)"     \
       "0.001*exp(-x^2)"  "0"               "0.001*exp(-x^2)"     \
       "0.005*exp(-x^2)"  "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"

QS_1 = "0.03*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))+0.03" "0.005*exp(-(x-7)^2)"                                "0.005*exp(-(x+7)^2)"                                   "0"                                              \
       "0.005*exp(-(x-7)^2)"                         "0.03*(tanh(1.6*(x-7))+tanh(1.6*x)+tanh(1.6*(x+7)))" "0.005*exp(-x^2)"                                       "0.005*exp(-(x+7)^2)"                            \
       "0.005*exp(-(x+7)^2)"                         "0.005*exp(-x^2)"                                    "(-0.03)*(tanh(1.6*x)+tanh(1.6*(x-7))+tanh(1.6*(x+7)))" "0.005*exp(-(x-7)^2)"                            \
       "0"                                           "0.005*exp(-(x+7)^2)"                                "0.005*exp(-(x-7)^2)"                                   "(-0.03)*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))-0.03" \

QS_2 = "0.01*tanh(0.6*x)" "0.001*exp(-x^2)"   "0.001*exp(-x^2)"      "0.005*exp(-x^2)"     \
       "0.001*exp(-x^2)"  "0.003*tanh(0.6*x)" "0.005*exp(-x^2)"      "0.001*exp(-x^2)"     \
       "0.001*exp(-x^2)"  "0.005*exp(-x^2)"   "(-0.003)*tanh(0.6*x)" "0.001*exp(-x^2)"     \
       "0.005*exp(-x^2)"  "0.001*exp(-x^2)"   "0.001*exp(-x^2)"      "(-0.01)*tanh(0.6*x)"

P = 15

random_cdyn_lz_1d_tully_1:
	"../bin/acorn_cdyn" -c -10 -e 0 -i 3500 -l 500 -m 2000 -p $(P) -r 1 -s 1 -t 1000 -u $(TULLY_1) --adiabatic

random_cdyn_lz_1d_tully_2:
	"../bin/acorn_cdyn" -c -10 -e 0 -i 3500 -l 500 -m 2000 -p $(P) -r 1 -s 1 -t 1000 -u $(TULLY_2) --adiabatic

random_cdyn_lz_1d_ds_1:
	"../bin/acorn_cdyn" -c -10 -e 0 -i 3500 -l 500 -m 2000 -p $(P) -r 1 -s 1 -t 1000 -u $(DS_1) --adiabatic

random_cdyn_lz_1d_ts_1:
	"../bin/acorn_cdyn" -c -10 -e 1 -i 3500 -l 500 -m 2000 -p $(P) -r 1 -s 1 -t 1000 -u $(TS_1) --adiabatic

random_cdyn_lz_1d_ts_2:
	"../bin/acorn_cdyn" -c -10 -e 0 -i 3500 -l 500 -m 2000 -p $(P) -r 1 -s 1 -t 1000 -u $(TS_2) --adiabatic

random_cdyn_lz_1d_qs_1:
	"../bin/acorn_cdyn" -c -10 -e 2 -i 3500 -l 500 -m 2000 -p $(P) -r 1 -s 1 -t 1000 -u $(QS_1) --adiabatic

random_cdyn_lz_1d_qs_2:
	"../bin/acorn_cdyn" -c -10 -e 0 -i 3500 -l 500 -m 2000 -p $(P) -r 1 -s 1 -t 1000 -u $(QS_2) --adiabatic

random_qdyn_1d_tully_1:
	"../bin/acorn_expression" -g -24 40 -o U_DIA.mat -p 1024 -e $(TULLY_1)
	"../bin/acorn_expression" -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	"../bin/acorn_qdyn" -i 350 -m 2000 -p $(P) -s 10 --adiabatic --savewfn --align --factor 0.01

random_qdyn_1d_tully_2:
	"../bin/acorn_expression" -g -24 40 -o U_DIA.mat -p 1024 -e $(TULLY_2)
	"../bin/acorn_expression" -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	"../bin/acorn_qdyn" -i 350 -m 2000 -p $(P) -s 10 --adiabatic --savewfn --align --factor 0.01

random_qdyn_1d_ds_1:
	"../bin/acorn_expression" -g -24 40 -o U_DIA.mat -p 1024 -e $(DS_1)
	"../bin/acorn_expression" -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0"
	"../bin/acorn_qdyn" -i 350 -m 2000 -p $(P) -s 10 --adiabatic --savewfn --align --factor 0.01

random_qdyn_1d_ts_1:
	"../bin/acorn_expression" -g -24 40 -o U_DIA.mat -p 1024 -e $(TS_1)
	"../bin/acorn_expression" -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "0" "0" "exp(-(x+10)^2)" "0" "0" "0"
	"../bin/acorn_qdyn" -i 350 -m 2000 -p $(P) -s 10 --adiabatic --savewfn --align --factor 0.05

random_qdyn_1d_ts_2:
	"../bin/acorn_expression" -g -24 40 -o U_DIA.mat -p 1024 -e $(TS_2)
	"../bin/acorn_expression" -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0" "0" "0"
	"../bin/acorn_qdyn" -i 350 -m 2000 -p $(P) -s 10 --adiabatic --savewfn --align --factor 0.01

random_qdyn_1d_qs_1:
	"../bin/acorn_expression" -g -24 40 -o U_DIA.mat -p 1024 -e $(QS_1)
	"../bin/acorn_expression" -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "0" "0" "0" "0" "exp(-(x+10)^2)" "0" "0" "0"
	"../bin/acorn_qdyn" -i 350 -m 2000 -p $(P) -s 10 --adiabatic --savewfn --align --factor 0.05

random_qdyn_1d_qs_2:
	"../bin/acorn_expression" -g -24 40 -o U_DIA.mat -p 1024 -e $(QS_2)
	"../bin/acorn_expression" -g -24 40 -o PSI_DIA_GUESS.mat -p 1024 -e "exp(-(x+10)^2)" "0" "0" "0" "0" "0" "0" "0"
	"../bin/acorn_qdyn" -i 350 -m 2000 -p $(P) -s 10 --adiabatic --savewfn --align --factor 0.01

random_qdyn_hydrogen_imaginary:
	"../bin/acorn_expression" -d 3 -g -16 16 -o U_DIA.mat -p 128 -e "(-1)/sqrt(x^2+y^2+z^2)"
	"../bin/acorn_expression" -d 3 -g -16 16 -o PSI_DIA_GUESS.mat -p 128 -e "exp(-(x^2+y^2+z^2))" "0"
	"../bin/acorn_qdyn" -d 3 -i 1000 -m 0.999455 -o 1 -p 0 -s 0.1

random_qdyn_kepler_real:
	"../bin/acorn_expression" -d 2 -g -16 16 -o U_DIA.mat -p 64 -e "10/(x^2+y^2)-10/sqrt(x^2+y^2)"
	"../bin/acorn_expression" -d 2 -g -16 16 -o PSI_DIA_GUESS.mat -p 64 -e "exp(-((x-1)^2+(y-1)^2))" "0"
	"../bin/acorn_qdyn" -d 2 -i 1000 -m 1 -p 0 -s 0.1 --savewfn

random_dyn_1d_tully_1: random_qdyn_1d_tully_1 random_cdyn_lz_1d_tully_1
random_dyn_1d_tully_2: random_qdyn_1d_tully_2 random_cdyn_lz_1d_tully_2
random_dyn_1d_ds_1: random_qdyn_1d_ds_1 random_cdyn_lz_1d_ds_1
random_dyn_1d_ts_1: random_qdyn_1d_ts_1 random_cdyn_lz_1d_ts_1
random_dyn_1d_ts_2: random_qdyn_1d_ts_2 random_cdyn_lz_1d_ts_2
random_dyn_1d_qs_1: random_qdyn_1d_qs_1 random_cdyn_lz_1d_qs_1
random_dyn_1d_qs_2: random_qdyn_1d_qs_2 random_cdyn_lz_1d_qs_2
