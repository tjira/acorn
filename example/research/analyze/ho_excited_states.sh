#!/bin/bash

$PLOT SPECTRUM.mat --xlabel "Energy (a.u)" --ylabel "Intensity" --output SPECTRUM.png --png
$PLOT ACF.mat      --xlabel "Time (a.u)"   --ylabel "Overlap"   --output      ACF.png --png

$PLOT WAVEFUNCTION.mat:0,1 --animate 2 --fps 300 --xlabel "Coordinate (a.u)" --ylabel "Value" --output WAVEFUNCTION --mp4
