#!/bin/bash

# run the dynamics
acorn_cdyn -d "x" -i 1000 -m 1 -p 0 -s 0.01 -t 1 -u "0.5*x^2"
