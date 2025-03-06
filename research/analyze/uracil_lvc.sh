#!/bin/bash

$PLOT POPULATION_FSSH_8D.mat:0,1,2 resource/uracil_lvc/8D_D0_AFSSH-94.dat resource/uracil_lvc/8D_D1_AFSSH-94.dat resource/uracil_lvc/8D_D2_AFSSH-94.dat --legend "D0 (FSSH)" "D1 (FSSH)" "D2 (FSSH)" "D0 (REF A-FSSH)" "D1 (REF A-FSSH)" "D2 (REF A-FSSH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_8D_REFERENCE_AFSSH.png

$PLOT POPULATION_FSSH_8D.mat:0 POPULATION_LZSH_8D.mat:0 POPULATION_KTSH_8D.mat:0 POPULATION_MASH_8D.mat:0 resource/uracil_lvc/8D_D0_MCTDH.dat --legend "D0 (FSSH)" "D0 (LZSH)" "D0 (KTSH)" "D0 (MASH)" "D0 (MCTDH)" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION_8D_D0_COMPARISON.png

$PLOT ADIABATIC_POTENTIAL_12D_Q3.mat:0,5,10,15  --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{13}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q3.png
$PLOT ADIABATIC_POTENTIAL_12D_Q7.mat:0,5,10,15  --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{17}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q7.png
$PLOT ADIABATIC_POTENTIAL_12D_Q10.mat:0,5,10    --legend "D0" "D1" "D2"      --xlabel "\$q_{10}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q10.png
$PLOT ADIABATIC_POTENTIAL_12D_Q11.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{11}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q11.png
$PLOT ADIABATIC_POTENTIAL_12D_Q12.mat:0,5,10    --legend "D0" "D1" "D2"      --xlabel "\$q_{12}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q12.png
$PLOT ADIABATIC_POTENTIAL_12D_Q18.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{18}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q18.png
$PLOT ADIABATIC_POTENTIAL_12D_Q19.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{19}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q19.png
$PLOT ADIABATIC_POTENTIAL_12D_Q20.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{20}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q20.png
$PLOT ADIABATIC_POTENTIAL_12D_Q21.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{21}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q21.png
$PLOT ADIABATIC_POTENTIAL_12D_Q24.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{24}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q24.png
$PLOT ADIABATIC_POTENTIAL_12D_Q25.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{25}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q25.png
$PLOT ADIABATIC_POTENTIAL_12D_Q26.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{26}$" --ylabel "Energy (a.u)" --output ADIABATIC_POTENTIAL_12D_Q26.png

montage ADIABATIC_POTENTIAL_12D_Q3.png ADIABATIC_POTENTIAL_12D_Q7.png ADIABATIC_POTENTIAL_12D_Q11.png ADIABATIC_POTENTIAL_12D_Q18.png ADIABATIC_POTENTIAL_12D_Q19.png ADIABATIC_POTENTIAL_12D_Q20.png ADIABATIC_POTENTIAL_12D_Q21.png ADIABATIC_POTENTIAL_12D_Q24.png ADIABATIC_POTENTIAL_12D_Q25.png ADIABATIC_POTENTIAL_12D_Q26.png ADIABATIC_POTENTIAL_12D_Q10.png ADIABATIC_POTENTIAL_12D_Q12.png -mode concatenate -tile 4x3 ADIABATIC_POTENTIAL_12D.png

rm ADIABATIC_POTENTIAL_12D_Q*.png
