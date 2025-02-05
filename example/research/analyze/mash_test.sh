#!/bin/bash

$PLOT POPULATION_EXACT.mat:1 POPULATION_FSSH.mat:1 POPULATION_MASH_FSSH.mat:1 --legend "S0 (EXACT)" "S0 (FSSH)" "S0 (MA-FSSH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_FSSH --dpi 100 --png
$PLOT POPULATION_EXACT.mat:1 POPULATION_MASH.mat:1 misc/mash_test/MASH_REFERENCE_TULLY2.dat --legend "S0 (EXACT)" "S0 (MASH)" "S0 (MASH REF)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH_REFERENCE --dpi 100 --png
$PLOT POPULATION_EXACT.mat POPULATION_MASH.mat --legend "S0 (EXACT)" "S1 (EXACT)" "S0 (MASH)" "S1 (MASH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH --dpi 100 --png
