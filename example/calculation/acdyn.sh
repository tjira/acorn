#!/bin/bash

# run the dynamics
acorn_cdyn -c 3 -i 1000 -m 1 -p 0 -r 1 -s 0.01 -t 5 \
           -u "0.5*x^2" \
           -d "x"
