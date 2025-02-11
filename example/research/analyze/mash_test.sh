#!/bin/bash

$PLOT POPULATION_EXACT_tully1D_2.mat:1 POPULATION_FSSH_tully1D_2.mat:1 POPULATION_MASH_FSSH_tully1D_2.mat:1 --legend "S0 (EXACT)" "S0 (FSSH)" "S0 (MA-FSSH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_FSSH_tully1D_2 --png
$PLOT POPULATION_EXACT_tully1D_2.mat:1 POPULATION_MASH_tully1D_2.mat:1 misc/mash_test/MASH_REFERENCE_TULLY2.dat --legend "S0 (EXACT)" "S0 (MASH)" "S0 (MASH REF)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH_REFERENCE --png
$PLOT POPULATION_EXACT_tully1D_2.mat POPULATION_MASH_tully1D_2.mat --legend "S0 (EXACT)" "S1 (EXACT)" "S0 (MASH)" "S1 (MASH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_MASH_tully1D_2 --png

$PLOT POPULATION_EXACT_tripleState1D_3.mat POPULATION_FSSH_tripleState1D_3.mat POPULATION_MASH_tripleState1D_3.mat --legend "S0 (EXACT)" "S1 (EXACT)" "S2 (EXACT)" "S0 (FSSH)" "S1 (FSSH)" "S2 (FSSH)" "S0 (MASH)" "S1 (MASH)" "S2 (MASH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_FSSH_tripleState1D_3 --png
