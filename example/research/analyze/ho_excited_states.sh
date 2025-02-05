#!/bin/bash

$PLOT SPECTRUM_1D.mat --domain 0  6 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_1D --png
$PLOT SPECTRUM_2D.mat --domain 0  8 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_2D --png
$PLOT SPECTRUM_3D.mat --domain 0 10 --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM_3D --png

$PLOT ACF_1D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_1D --png
$PLOT ACF_2D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_2D --png
$PLOT ACF_3D.mat --xlabel "Time (a.u)" --ylabel "Autocorrelation Function" --output ACF_3D --png

# $PLOT WAVEFUNCTION_1D_REAL.mat:0,1 --animate 2 --fps 300 --xlabel "Coordinate (a.u)" --ylabel "Value" --output WAVEFUNCTION_1D --mp4
