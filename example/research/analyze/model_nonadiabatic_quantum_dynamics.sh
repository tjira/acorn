#!/bin/bash

$PLOT WAVEFUNCTION_DIABATIC_tully1D_1.mat:0,1,2,3~POTENTIAL_DIABATIC_tully1D_1.mat:0-0,1-0,2-3,3-3   --animate 4 --scale 0.01 --xlabel "Coordinate (a.u)" --ylabel "Value" --output WAVEFUNCTION_DIABATIC_tully1D_1  --gif
$PLOT WAVEFUNCTION_ADIABATIC_tully1D_1.mat:0,1,2,3~POTENTIAL_ADIABATIC_tully1D_1.mat:0-0,1-0,2-3,3-3 --animate 4 --scale 0.01 --xlabel "Coordinate (a.u)" --ylabel "Value" --output WAVEFUNCTION_ADIABATIC_tully1D_1 --gif

$PLOT WAVEFUNCTION_DIABATIC_tripleState1D_3.mat:0,1,2,3,4,5~POTENTIAL_DIABATIC_tripleState1D_3.mat:0-0,1-0,2-4,3-4,4-8,5-8   --animate 6 --scale 0.01 --xlabel "Coordinate (a.u)" --ylabel "Value" --output WAVEFUNCTION_DIABATIC_tripleState1D_3  --gif
$PLOT WAVEFUNCTION_ADIABATIC_tripleState1D_3.mat:0,1,2,3,4,5~POTENTIAL_ADIABATIC_tripleState1D_3.mat:0-0,1-0,2-4,3-4,4-8,5-8 --animate 6 --scale 0.01 --xlabel "Coordinate (a.u)" --ylabel "Value" --output WAVEFUNCTION_ADIABATIC_tripleState1D_3 --gif

$PLOT POPULATION_DIABATIC_tully1D_1.mat  --xlabel "Time (a.u)" --ylabel "Population" --legend "S0" "S1" --output POPULATION_DIABATIC_tully1D_1  --png
$PLOT POPULATION_ADIABATIC_tully1D_1.mat --xlabel "Time (a.u)" --ylabel "Population" --legend "S0" "S1" --output POPULATION_ADIABATIC_tully1D_1 --png

$PLOT POPULATION_DIABATIC_tripleState1D_3.mat  --xlabel "Time (a.u)" --ylabel "Population" --legend "S0" "S1" "S2" --output POPULATION_DIABATIC_tripleState1D_3  --png
$PLOT POPULATION_ADIABATIC_tripleState1D_3.mat --xlabel "Time (a.u)" --ylabel "Population" --legend "S0" "S1" "S2" --output POPULATION_ADIABATIC_tripleState1D_3 --png

montage POPULATION_DIABATIC_tully1D_1.png       POPULATION_ADIABATIC_tully1D_1.png       -mode concatenate -tile x1 POPULATION_tully1D_1.png
montage POPULATION_DIABATIC_tripleState1D_3.png POPULATION_ADIABATIC_tripleState1D_3.png -mode concatenate -tile x1 POPULATION_tripleState1D_3.png

ffmpeg -i WAVEFUNCTION_DIABATIC_tripleState1D_3.gif -i WAVEFUNCTION_ADIABATIC_tripleState1D_3.gif -filter_complex "palettegen" -loglevel error palette.png && ffmpeg -i WAVEFUNCTION_DIABATIC_tripleState1D_3.gif -i WAVEFUNCTION_ADIABATIC_tully1D_1.gif -i palette.png -filter_complex "[0:v][1:v]hstack[v]; [v][2:v]paletteuse" -loglevel error WAVEFUNCTION_tully1D_1.gif && rm palette.png

ffmpeg -i WAVEFUNCTION_DIABATIC_tripleState1D_3.gif -i WAVEFUNCTION_ADIABATIC_tripleState1D_3.gif -filter_complex "palettegen" -loglevel error palette.png && ffmpeg -i WAVEFUNCTION_DIABATIC_tripleState1D_3.gif -i WAVEFUNCTION_ADIABATIC_tripleState1D_3.gif -i palette.png -filter_complex "[0:v][1:v]hstack[v]; [v][2:v]paletteuse" -loglevel error WAVEFUNCTION_tripleState1D_3.gif && rm palette.png

rm *DIABATIC*.gif *DIABATIC*.png
