#!/bin/bash

$PLOT POPULATION_EXACT_tully1D_2.mat:1 POPULATION_MASH_tully1D_2.mat:1 misc/mash_test/MASH_REFERENCE_TULLY2.dat --legend "S0 (EXACT)" "S0 (MASH)" "S0 (MASH REF)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH_REFERENCE_tully1D_2 --png
$PLOT POPULATION_EXACT_tully1D_2.mat POPULATION_MASH_tully1D_2.mat --legend "S0 (EXACT)" "S1 (EXACT)" "S0 (MASH)" "S1 (MASH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH_tully1D_2 --png

$PLOT POPULATION_EXACT_tully1D_3.mat POPULATION_MASH_tully1D_3.mat --legend "S0 (EXACT)" "S1 (EXACT)" "S0 (MASH)" "S1 (MASH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH_tully1D_3 --png
