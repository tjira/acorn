#!/bin/bash

$PLOT ACF_1D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_1D.png
$PLOT ACF_2D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_2D.png

$PLOT ACF_TRANSFORMED_1D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_TRANSFORMED_1D.png
$PLOT ACF_TRANSFORMED_2D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_TRANSFORMED_2D.png

$PLOT SPECTRUM_1D.mat --xlim 0 6 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_1D.png
$PLOT SPECTRUM_2D.mat --xlim 0 8 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_2D.png

$PLOT SPECTRUM_TRANSFORMED_1D.mat --xlim 0 6 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_TRANSFORMED_1D.png
$PLOT SPECTRUM_TRANSFORMED_2D.mat --xlim 0 8 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_TRANSFORMED_2D.png

$PLOT WAVEFUNCTION_1D_REAL.mat:0,1 --xlabel "Coordinate (a.u)" --ylabel "Value" --animate 2 --output WAVEFUNCTION_1D.gif

lines.py POTENTIAL_1D.mat WAVEFUNCTION_1D_IMAG_1.mat:2000,2001 WAVEFUNCTION_1D_IMAG_2.mat:2000,2001 WAVEFUNCTION_1D_IMAG_3.mat:2000,2001 WAVEFUNCTION_1D_IMAG_4.mat:2000,2001 WAVEFUNCTION_1D_IMAG_5.mat:2000,2001 WAVEFUNCTION_1D_IMAG_6.mat:2000,2001 WAVEFUNCTION_1D_IMAG_7.mat:2000,2001 WAVEFUNCTION_1D_IMAG_8.mat:2000,2001 WAVEFUNCTION_1D_IMAG_9.mat:2000,2001 --xlabel "Coordinate (a.u)" --ylabel "Energy (a.u)" --xlim -5 5 --ylim 0 9 --offsets 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 0.5 0.5 1.5 1.5 2.5 2.5 3.5 3.5 4.5 4.5 5.5 5.5 6.5 6.5 7.5 7.5 8.5 8.5 --scales 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 -o EIGENSTATES_1D.png

montage ACF_1D.png SPECTRUM_1D.png -mode concatenate -tile x1 JOIN_ACF_SPECTRUM_1D.png
montage ACF_2D.png SPECTRUM_2D.png -mode concatenate -tile x1 JOIN_ACF_SPECTRUM_2D.png

montage ACF_TRANSFORMED_1D.png SPECTRUM_TRANSFORMED_1D.png -mode concatenate -tile x1 JOIN_ACF_SPECTRUM_TRANSFORMED_1D.png
montage ACF_TRANSFORMED_2D.png SPECTRUM_TRANSFORMED_2D.png -mode concatenate -tile x1 JOIN_ACF_SPECTRUM_TRANSFORMED_2D.png

rm ACF_*.png SPECTRUM_*.png
