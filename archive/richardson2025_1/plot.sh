#!/bin/bash

lines.py POTENTIAL.mat:0,3 POPULATION_EXACT.mat:1 POPULATION_MEAN_FSSH.mat:1 POPULATION_MEAN_LZSH.mat:1 POPULATION_MEAN_MASH.mat:1 --figsize 4 12 --legends every "S\$_0\$" "S\$_1\$" "EXACT" "FSSH" "LZSH" "MASH" --subplots 121 121 122 122 122 122 122 --xlabel "Coordinate (a.u.)" "Time (a.u.)" --ylabel "Energy (a.u.)" "S\$_1\$ Populations" --ylim nan nan -0.02 1.02 --output POPULATIONS.png

lines.py BLOCH_VECTOR_MEAN_FSSH.mat BLOCH_VECTOR_MEAN_MASH.mat BLOCH_VECTOR_EXACT.mat --figsize 4 18 --legends every "\$<\sigma_x>\$" "\$<\sigma_y>\$" "\$<\sigma_z>\$" "\$<\sigma_x>\$" "\$<\sigma_y>\$" "\$<\sigma_z>\$" "\$<\sigma_x>\$" "\$<\sigma_y>\$" "\$<\sigma_z>\$" --subplots 131 131 131 132 132 132 133 133 133 --title "FSSH" "MASH" "EXACT" --xlabel "Time(a.u.)" "Time(a.u.)" "Time(a.u.)" --ylabel "Bloch Vector Element" "Bloch Vector Element" "Bloch Vector Element" --ylim -1.02 1.02 -1.02 1.02 -1.02 1.02 --output S.png
