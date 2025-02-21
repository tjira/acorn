#!/bin/bash

$PLOT ACF_1D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_1D --png
$PLOT ACF_2D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_2D --png

$PLOT ACF_TRANSFORMED_1D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_TRANSFORMED_1D --png
$PLOT ACF_TRANSFORMED_2D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_TRANSFORMED_2D --png

$PLOT SPECTRUM_1D.mat --domain 0  6 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_1D --png
$PLOT SPECTRUM_2D.mat --domain 0  8 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_2D --png

$PLOT SPECTRUM_TRANSFORMED_1D.mat --domain 0  6 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_TRANSFORMED_1D --png
$PLOT SPECTRUM_TRANSFORMED_2D.mat --domain 0  8 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_TRANSFORMED_2D --png

$PLOT WAVEFUNCTION_1D_REAL.mat:0,1 --animate 2 --ndim 1 --xlabel "Coordinate (a.u)"    --ylabel "Value"               --output WAVEFUNCTION_1D --gif
$PLOT WAVEFUNCTION_2D_REAL.mat:0,1 --animate 2 --ndim 2 --xlabel "Coordinate #1 (a.u)" --ylabel "Coordinate #2 (a.u)" --output WAVEFUNCTION_2D --gif

$PLOT POTENTIAL_1D.mat WAVEFUNCTION_1D_IMAG_1.mat:1000,1001 WAVEFUNCTION_1D_IMAG_2.mat:1000,1001 WAVEFUNCTION_1D_IMAG_3.mat:1000,1001 WAVEFUNCTION_1D_IMAG_4.mat:1000,1001 WAVEFUNCTION_1D_IMAG_5.mat:1000,1001 WAVEFUNCTION_1D_IMAG_6.mat:1000,1001 WAVEFUNCTION_1D_IMAG_7.mat:1000,1001 WAVEFUNCTION_1D_IMAG_8.mat:1000,1001 WAVEFUNCTION_1D_IMAG_9.mat:1000,1001 --offset 0 0.5 0.5 1.5 1.5 2.5 2.5 3.5 3.5 4.5 4.5 5.5 5.5 6.5 6.5 7.5 7.5 8.5 8.5 --xlabel "Coordinate (a.u)" --ylabel "Energy (a.u)" --domain -5 5 --value 0 9 --scale 0.8 --output EIGENSTATES_1D --png

$PLOT WAVEFUNCTION_2D_IMAG_1.mat:1000,1001 WAVEFUNCTION_2D_IMAG_2.mat:1000,1001 WAVEFUNCTION_2D_IMAG_3.mat:1000,1001 WAVEFUNCTION_2D_IMAG_4.mat:1000,1001 WAVEFUNCTION_2D_IMAG_5.mat:1000,1001 WAVEFUNCTION_2D_IMAG_6.mat:1000,1001 WAVEFUNCTION_2D_IMAG_7.mat:1000,1001 WAVEFUNCTION_2D_IMAG_8.mat:1000,1001 WAVEFUNCTION_2D_IMAG_9.mat:1000,1001 --ndim 2 --blank --output EIGENSTATES_2D --png

montage ACF_1D.png SPECTRUM_1D.png -mode concatenate -tile x1 JOIN_ACF_SPECTRUM_1D.png
montage ACF_2D.png SPECTRUM_2D.png -mode concatenate -tile x1 JOIN_ACF_SPECTRUM_2D.png

montage ACF_TRANSFORMED_1D.png SPECTRUM_TRANSFORMED_1D.png -mode concatenate -tile x1 JOIN_ACF_SPECTRUM_TRANSFORMED_1D.png
montage ACF_TRANSFORMED_2D.png SPECTRUM_TRANSFORMED_2D.png -mode concatenate -tile x1 JOIN_ACF_SPECTRUM_TRANSFORMED_2D.png

rm ACF_*.png SPECTRUM_*.png
