#!/bin/bash

# generate the potential matrix
acorn_expression -e "0.5*x^2" -g -16 16 -p 2048 -o U_ADIA.mat

# generate the initial wavefunction guess
acorn_expression -e "exp(-(x-1)^2)" "0" -g -16 16 -p 2048 -o PSI_ADIA_GUESS.mat

# run the dynamics
acorn_aqdyn -i 200 -m 1 -n 3 -p 0 -s 0.1 --imaginary
