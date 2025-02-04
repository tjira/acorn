#!/bin/bash

$PLOT POPULATION_EXACT.mat:1 POPULATION_FSSH.mat:1 POPULATION_MASH.mat:1 --legend "S0 (EXACT)" "S0 (FSSH)" "S0 (MASH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION --dpi 100 --png
