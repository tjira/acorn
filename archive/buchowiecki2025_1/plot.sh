#!/bin/bash

lines.py POTENTIAL_DIA_1D.mat:0,1,3 POTENTIAL_ADIA_1D.mat:0,3 --figsize 3 12 --subplots 131 132 131 133 133 --title "Diabatic Potential" "SOC" "Adiabatic Potential" --xlabel "Coordinate (a.u.)" "Coordinate (a.u.)" "Coordinate (a.u.)" --xlim 1 10 0.001 15 1 10 --ylabel "Energy (a.u.)" "SOC (a.u.)" "Energy (a.u.)" --ylim -0.2 0.2 nan nan -0.2 0.2 --output POTENTIAL_1D.png

# echo $(ls -l POPULATION_H_1D_E=*.mat | wc -l)" 4" > PROBABILITY_EXACT_1D.mat
#
# for POP_H in $(ls POPULATION_H_1D_E=*.mat); do
#
#     E=$(echo "$POP_H" | sed 's/POPULATION_H_1D_E=// ; s/.mat//')
#
#     POP_D=$(echo "$POP_H" | sed 's/H/D/')
#     POP_T=$(echo "$POP_H" | sed 's/H/T/')
#
#     TIME_H=$(cat $(echo "$POP_H" | sed 's/POPULATION/POSITION/') | awk '$1 > 0 && $2 > 15 { print $1; exit }')
#     TIME_D=$(cat $(echo "$POP_D" | sed 's/POPULATION/POSITION/') | awk '$1 > 0 && $2 > 15 { print $1; exit }')
#     TIME_T=$(cat $(echo "$POP_T" | sed 's/POPULATION/POSITION/') | awk '$1 > 0 && $2 > 15 { print $1; exit }')
#
#     PROB_H=$(cat $POP_H | awk -v TIME=$TIME_H '$1 == TIME { print $3; exit }')
#     PROB_D=$(cat $POP_D | awk -v TIME=$TIME_D '$1 == TIME { print $3; exit }')
#     PROB_T=$(cat $POP_T | awk -v TIME=$TIME_T '$1 == TIME { print $3; exit }')
#
#     echo "$E $PROB_H $PROB_D $PROB_T" >> PROBABILITY_EXACT_1D.mat
# done
#
# matsort -i PROBABILITY_EXACT_1D.mat -o PROBABILITY_EXACT_1D.mat

lines.py P_LZ.dat PROBABILITY_EXACT_1D.mat P_LZ.dat PROBABILITY_EXACT_1D.mat P_LZ.dat PROBABILITY_EXACT_1D.mat P_LZ.dat PROBABILITY_EXACT_1D.mat --alphas 0,1,2,6,7,8,12,13,14,18,19,20 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 --colors 3,4,5,9,10,11,15,16,17,21,22,23 "darkblue" "darkorange" "darkgreen" "darkblue" "darkorange" "darkgreen" "darkblue" "darkorange" "darkgreen" "darkblue" "darkorange" "darkgreen" "darkblue" "darkorange" "darkgreen" --figsize 6 8 --legends every "H (LZ)" "D (LZ)" "T (LZ)" "H (QD)" "D (QD)" "T (QD)" "H (LZ)" "D (LZ)" "T (LZ)" "H (QD)" "D (QD)" "T (QD)" "H (LZ)" "D (LZ)" "T (LZ)" "H (QD)" "D (QD)" "T (QD)" "H (LZ)" "D (LZ)" "T (LZ)" "H (QD)" "D (QD)" "T (QD)" --legpos 0.6 1 0.6 1 0.6 1 0.6 1 --styles 3,4,5,9,10,11,15,16,17,21,22,23,27,28,29 dashed dashed dashed dashed dashed dashed dashed dashed dashed dashed dashed dashed --subplots 221 221 221 221 221 221 222 222 222 222 222 222 223 223 223 223 223 223 224 224 224 224 224 224 --xlabel "Kinetic Energy (a.u.)" "Kinetic Energy (a.u.)" "Kinetic Energy (a.u.)" --xlim 0.0347289 0.1 0.1 1 1 10 10 100 --ylabel "Charge Transfer Probability" "Charge Transfer Probability" "Charge Transfer Probability" --ylim 0 0.035 0 0.01 0 0.0035 0 0.0007 --output PROBABILITY_SPLIT_1D.png

lines.py P_LZ.dat PROBABILITY_EXACT_1D.mat --alphas 0,1,2 0.7 0.7 0.7 --colors 3,4,5 "darkblue" "darkorange" "darkgreen" --legends every "H (LZ)" "D (LZ)" "T (LZ)" "H (QD)" "D (QD)" "T (QD)" --legpos 0.3 1 --styles 3,4,5 dashed dashed dashed --xlabel "Kinetic Energy (a.u.)" --xlim 0.0347289 100 --ylabel "Charge Transfer Probability" --ylim 0 0.001 --output PROBABILITY_FULL_1D.png
