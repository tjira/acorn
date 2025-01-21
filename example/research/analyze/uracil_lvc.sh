#!/bin/bash

$PLOT POPULATION_FSSH.mat --legend "D0" "D1" "D2" "D3" --xlabel "Time (a.u.)" --ylabel "Population" --output POPULATION --png

$PLOT POTENTIAL_Q10.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{10}$" --ylabel "Energy (a.u)" --output POTENTIAL_Q10 --png --dpi 100
$PLOT POTENTIAL_Q12.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{12}$" --ylabel "Energy (a.u)" --output POTENTIAL_Q12 --png --dpi 100
$PLOT POTENTIAL_Q18.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{18}$" --ylabel "Energy (a.u)" --output POTENTIAL_Q18 --png --dpi 100
$PLOT POTENTIAL_Q20.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{20}$" --ylabel "Energy (a.u)" --output POTENTIAL_Q20 --png --dpi 100
$PLOT POTENTIAL_Q21.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{21}$" --ylabel "Energy (a.u)" --output POTENTIAL_Q21 --png --dpi 100
$PLOT POTENTIAL_Q24.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{24}$" --ylabel "Energy (a.u)" --output POTENTIAL_Q24 --png --dpi 100
$PLOT POTENTIAL_Q25.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{25}$" --ylabel "Energy (a.u)" --output POTENTIAL_Q25 --png --dpi 100
$PLOT POTENTIAL_Q26.mat:0,5,10,15 --legend "D0" "D1" "D2" "D3" --xlabel "\$q_{26}$" --ylabel "Energy (a.u)" --output POTENTIAL_Q26 --png --dpi 101

montage POTENTIAL_Q10.png POTENTIAL_Q12.png POTENTIAL_Q18.png POTENTIAL_Q20.png POTENTIAL_Q21.png POTENTIAL_Q24.png POTENTIAL_Q25.png POTENTIAL_Q26.png -mode concatenate -tile 4x2 POTENTIAL.png
