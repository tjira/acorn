#!/bin/bash

$PLOT POPULATION_EXACT.mat:1 POPULATION_MASH.mat:1 misc/mash_test/MASH_REFERENCE_TULLY2.dat --legend "S0 (EXACT)" "S0 (MASH)" "S0 (MASH REF)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH_REFERENCE --png
$PLOT POPULATION_EXACT.mat POPULATION_MASH.mat --legend "S0 (EXACT)" "S1 (EXACT)" "S0 (MASH)" "S1 (MASH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH --png
