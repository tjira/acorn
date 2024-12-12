#!/bin/bash

PLOT_1D=./../../script/plot-1d.py

$PLOT_1D POTENTIAL_Q10.mat:0,5,10,15 --xlabel "\$q_{10}$" --ylabel "Energy (a.u.)" --legend "\$D_0$" "\$D_1$" "\$D_2$" "\$D_3$" --output POTENTIAL_Q10 --dpi 100 --png
$PLOT_1D POTENTIAL_Q12.mat:0,5,10,15 --xlabel "\$q_{12}$" --ylabel "Energy (a.u.)" --legend "\$D_0$" "\$D_1$" "\$D_2$" "\$D_3$" --output POTENTIAL_Q12 --dpi 100 --png
$PLOT_1D POTENTIAL_Q18.mat:0,5,10,15 --xlabel "\$q_{18}$" --ylabel "Energy (a.u.)" --legend "\$D_0$" "\$D_1$" "\$D_2$" "\$D_3$" --output POTENTIAL_Q18 --dpi 100 --png
$PLOT_1D POTENTIAL_Q20.mat:0,5,10,15 --xlabel "\$q_{20}$" --ylabel "Energy (a.u.)" --legend "\$D_0$" "\$D_1$" "\$D_2$" "\$D_3$" --output POTENTIAL_Q20 --dpi 100 --png
$PLOT_1D POTENTIAL_Q21.mat:0,5,10,15 --xlabel "\$q_{21}$" --ylabel "Energy (a.u.)" --legend "\$D_0$" "\$D_1$" "\$D_2$" "\$D_3$" --output POTENTIAL_Q21 --dpi 100 --png
$PLOT_1D POTENTIAL_Q24.mat:0,5,10,15 --xlabel "\$q_{24}$" --ylabel "Energy (a.u.)" --legend "\$D_0$" "\$D_1$" "\$D_2$" "\$D_3$" --output POTENTIAL_Q24 --dpi 100 --png
$PLOT_1D POTENTIAL_Q25.mat:0,5,10,15 --xlabel "\$q_{25}$" --ylabel "Energy (a.u.)" --legend "\$D_0$" "\$D_1$" "\$D_2$" "\$D_3$" --output POTENTIAL_Q25 --dpi 100 --png
$PLOT_1D POTENTIAL_Q26.mat:0,5,10,15 --xlabel "\$q_{26}$" --ylabel "Energy (a.u.)" --legend "\$D_0$" "\$D_1$" "\$D_2$" "\$D_3$" --output POTENTIAL_Q26 --dpi 100 --png

montage POTENTIAL_Q10.png POTENTIAL_Q12.png POTENTIAL_Q18.png POTENTIAL_Q20.png POTENTIAL_Q21.png POTENTIAL_Q24.png POTENTIAL_Q25.png POTENTIAL_Q26.png -mode concatenate -tile 4x2 POTENTIAL.png
