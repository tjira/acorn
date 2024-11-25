#!/bin/bash

PLOT_1D=./../../script/plot/plot-1d.py; MODELS=($(ls -l exact* | grep -oP "(?<=exact_)[^_]+_[^_]+" | sort | uniq)); MOMENTA=($(ls -l exact* | grep -oP "P=\K[0-9.]+" | sed 's/.$//' | sort | uniq))

for MODEL in ${MODELS[@]}; do

    NSTATE=$(tail -n 1 POPULATION_${MODEL}_* | tail -n 1 | awk '{print NF - 1}'); LEGEND_POP=(); LEGEND_POT=(); INDICES_POT="";

    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} EXACT"); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} FSSH" ); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} KFSSH"); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} LZSH" ); done

    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POT+=("S${i}"); done; for (( i=0; i<$NSTATE; i++ )); do INDICES_POT+=$(echo "$i*$NSTATE+$i" | bc -l)","; done

    $PLOT_1D  "POTENTIAL_DIABATIC_${MODEL}.mat:${INDICES_POT::-1}" --title  "DIABATIC POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --legend "${LEGEND_POT[@]}" --output  "POTENTIAL_DIABATIC_${MODEL}" --png
    $PLOT_1D "POTENTIAL_ADIABATIC_${MODEL}.mat:${INDICES_POT::-1}" --title "ADIABATIC POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --legend "${LEGEND_POT[@]}" --output "POTENTIAL_ADIABATIC_${MODEL}" --png

    echo "${#MOMENTA[@]} 4" > "FINAL_POPULATION_P_${MODEL}.mat"

    for MOMENTUM in ${MOMENTA[@]}; do

        IS_COL_EXACT=$(head -n 2 "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        IS_COL_FSSH=$( head -n 2  "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        IS_COL_KFSSH=$(head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        IS_COL_LZSH=$( head -n 2  "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')

        POP_EXACT=$(tail -n 1 "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      | awk -v i=$IS_COL_EXACT '{print $(i+2)}');
        POP_FSSH=$( tail -n 1  "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" | awk -v i=$IS_COL_FSSH  '{print $(i+2)}');
        POP_KFSSH=$(tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" | awk -v i=$IS_COL_KFSSH '{print $(i+2)}');
        POP_LZSH=$( tail -n 1  "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" | awk -v i=$IS_COL_LZSH  '{print $(i+2)}');

        echo "${MOMENTUM} ${POP_EXACT} ${POP_FSSH} ${POP_KFSSH} ${POP_LZSH}" >> "FINAL_POPULATION_P_${MODEL}.mat"

        $PLOT_1D "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat" "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" --title "ADIABATIC POPULATION: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "State Population" --legend "${LEGEND_POP[@]}" --output "POPULATION_${MODEL}_P=${MOMENTUM}" --png
        $PLOT_1D "KINETIC_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat" "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" --title "KINETIC ENERGY: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "Energy (a.u.)" --legend "EXACT" "FSSH" "KFSSH" "LZSH" --output "KINETIC_ENERGY_${MODEL}_P=${MOMENTUM}" --png
        $PLOT_1D "POTENTIAL_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat" "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" --title "POTENTIAL ENERGY: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "Energy (a.u.)" --legend "EXACT" "FSSH" "KFSSH" "LZSH" --output "POTENTIAL_ENERGY_${MODEL}_P=${MOMENTUM}" --png
        $PLOT_1D "TOTAL_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat" "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" --title "TOTAL ENERGY: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "Energy (a.u.)" --legend "EXACT" "FSSH" "KFSSH" "LZSH" --output "TOTAL_ENERGY_${MODEL}_P=${MOMENTUM}" --png
        $PLOT_1D "POSITION_${MODEL}_P=${MOMENTUM}_EXACT.mat" "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" --title "POSITION: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "Coordinate (a.u.)" --legend "EXACT" "FSSH" "KFSSH" "LZSH" --output "POSITION_${MODEL}_P=${MOMENTUM}" --png
        $PLOT_1D "MOMENTUM_${MODEL}_P=${MOMENTUM}_EXACT.mat" "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" --title "MOMENTUM: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "Momentum (a.u.)" --legend "EXACT" "FSSH" "KFSSH" "LZSH" --output "MOMENTUM_${MODEL}_P=${MOMENTUM}" --png

        montage "POTENTIAL_ADIABATIC_${MODEL}.png" "POPULATION_${MODEL}_P=${MOMENTUM}.png" "POSITION_${MODEL}_P=${MOMENTUM}.png" "MOMENTUM_${MODEL}_P=${MOMENTUM}.png" -mode concatenate -tile x1 "TRAJECTORY_GENERAL_${MODEL}_P=${MOMENTUM}.png"
        montage "POTENTIAL_ADIABATIC_${MODEL}.png" "KINETIC_ENERGY_${MODEL}_P=${MOMENTUM}.png" "POTENTIAL_ENERGY_${MODEL}_P=${MOMENTUM}.png" "TOTAL_ENERGY_${MODEL}_P=${MOMENTUM}.png" -mode concatenate -tile x1 "TRAJECTORY_ENERGY_${MODEL}_P=${MOMENTUM}.png"
    done

    $PLOT_1D "FINAL_POPULATION_P_${MODEL}.mat" --title "MOMENTUM DEPENDENT IS POPULATION: ${MODEL}" --xlabel "Momentum (a.u.)" --ylabel "Initial State Population" --legend "EXACT" "FSSH" "KFSSH" "LZSH" --output "FINAL_POPULATION_P_${MODEL}" --png

    montage "POTENTIAL_DIABATIC_${MODEL}.png"  "POTENTIAL_ADIABATIC_${MODEL}.png" -mode concatenate -tile x1        "POTENTIAL_${MODEL}.png"
    montage "POTENTIAL_ADIABATIC_${MODEL}.png" "FINAL_POPULATION_P_${MODEL}.png"  -mode concatenate -tile x1 "FINAL_POPULATION_${MODEL}.png"
done

FILES=(${MODELS[@]/#/POTENTIAL_});        montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x        "AA_POTENTIAL.png"
FILES=(${MODELS[@]/#/FINAL_POPULATION_}); montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x "AA_FINAL_POPULATION.png"

for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORY_GENERAL_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "AA_TRAJECTORY_GENERAL_P=${MOMENTUM}.png"
done

for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORY_ENERGY_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "AA_TRAJECTORY_ENERGY_P=${MOMENTUM}.png"
done
