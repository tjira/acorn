#!/bin/bash

$PLOT POPULATION_EXACT.mat POPULATION_FSSH.mat POPULATION_MASH.mat --legend "S0 (EXACT)" "S1 (EXACT)" "S0 (FSSH)" "S1 (FSSH)" "S0 (MASH)" "S1 (MASH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION --png
