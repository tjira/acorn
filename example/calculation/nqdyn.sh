#!/bin/bash

# generate the potential matrix
acorn_expression -g -24 24 -p 4096 -o U_DIA.mat \
                 -e "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)"

# generate the initial wavefunction guess
acorn_expression -g -24 24 -p 4096 -o PSI_DIA_GUESS.mat \
                 -e "exp(-(x+10)^2)" "0" "0" "0"

# run the dynamics
acorn_nqdyn -i 350 -m 2000 -p 10.95 -s 10 --adiabatic
