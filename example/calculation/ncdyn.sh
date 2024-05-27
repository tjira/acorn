#!/bin/bash

# run the dynamics
acorn_cdyn -c -10 -i 350 -m 2000 -p 10.95 -r 1 -s 10 -t 100 \
           -u "0.01*tanh(0.6*x)" "0.001*exp(-x^2)" "0.001*exp(-x^2)" "(-0.01)*tanh(0.6*x)" \
           -d "0.006/cosh(0.6*x)^2" "(-0.002)*x*exp(-x^2)" "(-0.002)*x*exp(-x^2)" "(-0.006)/cosh(0.6*x)^2"
