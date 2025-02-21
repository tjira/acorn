#!/bin/bash

MODELS=($(ls -l POTENTIAL_DIABATIC* | grep -oP "(?<=POTENTIAL_DIABATIC_)[^_]+_[^_]" | sort | uniq)); MOMENTA=($(ls -l POPULATION* | grep -oP "P=\K[0-9.]+" | sort | uniq))

for MODEL in ${MODELS[@]}; do

    NSTATE=$(tail -n 1 POPULATION_${MODEL}_* | tail -n 1 | awk '{print NF - 1}');

    LEGEND_POP=(); LEGEND_FSSH=(); LEGEND_FSSH_COEF=(); LEGEND_KTSH_COEF=(); LEGEND_POT=(); LEGEND_TDC=(); INDICES_POT=""; INDICES_TDC="";

    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} EXACT"); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} FSSH" ); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} LZSH" ); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} MASH" ); done

    for (( i=0; i<$NSTATE; i++ )); do LEGEND_FSSH+=("S${i} FSSH"); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_FSSH+=("S${i} KTSH"); done

    for (( i=0; i<$NSTATE; i++ )); do LEGEND_FSSH_COEF+=("S${i} FSSH"        ); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_FSSH_COEF+=("S${i} FSSH (COEF)" ); done

    for (( i=0; i<$NSTATE; i++ )); do LEGEND_KTSH_COEF+=("S${i} KTSH"        ); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_KTSH_COEF+=("S${i} KTSH (COEF)" ); done

    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POT+=("S${i}"); done;

    for (( i=0; i<$NSTATE; i++ )); do for (( j=((i+1)); j<$NSTATE; j++ )); do LEGEND_TDC+=("T\$_{$i,$j}\$ (NUMERIC)"); done; done;
    for (( i=0; i<$NSTATE; i++ )); do for (( j=((i+1)); j<$NSTATE; j++ )); do LEGEND_TDC+=("T\$_{$i,$j}\$ (BAECKAN)"); done; done;

    for (( i=0; i<$NSTATE; i++ )); do INDICES_POT+=$(echo "$i*$NSTATE+$i" | bc -l)","; done

    for (( i=0; i<$NSTATE; i++ )); do for (( j=((i+1)); j<$NSTATE; j++ )); do INDICES_TDC+=$(echo "$i*$NSTATE+$j" | bc -l)","; done; done;

    $PLOT  "POTENTIAL_DIABATIC_${MODEL}.mat:${INDICES_POT::-1}" --title  "DIABATIC POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --legend "${LEGEND_POT[@]}" --output  "POTENTIAL_DIABATIC_${MODEL}" --png
    $PLOT "POTENTIAL_ADIABATIC_${MODEL}.mat:${INDICES_POT::-1}" --title "ADIABATIC POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --legend "${LEGEND_POT[@]}" --output "POTENTIAL_ADIABATIC_${MODEL}" --png

    echo "${#MOMENTA[@]} 4" > "FINAL_POPULATION_P_${MODEL}.mat"

    for MOMENTUM in ${MOMENTA[@]}; do

        [[ -f "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      ]] && IS_COL_EXACT=$(head -n 2 "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  ]] && IS_COL_FSSH=$( head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"  ]] && IS_COL_KTSH=$( head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat" | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  ]] && IS_COL_LZSH=$( head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"  ]] && IS_COL_MASH=$( head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"  | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')

        [[ -f "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      ]] && POP_EXACT=$(tail -n 1 "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      | awk -v i=$IS_COL_EXACT '{print $(i+2)}');
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  ]] && POP_FSSH=$( tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  | awk -v i=$IS_COL_FSSH  '{print $(i+2)}');
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"  ]] && POP_KTSH=$( tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"  | awk -v i=$IS_COL_KTSH  '{print $(i+2)}');
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  ]] && POP_LZSH=$( tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  | awk -v i=$IS_COL_LZSH  '{print $(i+2)}');
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"  ]] && POP_MASH=$( tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"  | awk -v i=$IS_COL_LZSH  '{print $(i+2)}');

        echo "${MOMENTUM} ${POP_EXACT} ${POP_FSSH} ${POP_KTSH} ${POP_LZSH}" >> "FINAL_POPULATION_P_${MODEL}.mat"

        $PLOT "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"       "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"       "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"       "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"       --title "ADIABATIC POPULATION: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "State Population"  --legend "${LEGEND_POP[@]}"           --output "POPULATION_${MODEL}_P=${MOMENTUM}"       --png
        $PLOT "KINETIC_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat"   "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"   "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"   "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"   --title "KINETIC ENERGY: ${MODEL}"       --xlabel "Time (a.u.)" --ylabel "Energy (a.u.)"     --legend "EXACT" "FSSH" "LZSH" "MASH" --output "KINETIC_ENERGY_${MODEL}_P=${MOMENTUM}"   --png
        $PLOT "POTENTIAL_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat" "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat" --title "POTENTIAL ENERGY: ${MODEL}"     --xlabel "Time (a.u.)" --ylabel "Energy (a.u.)"     --legend "EXACT" "FSSH" "LZSH" "MASH" --output "POTENTIAL_ENERGY_${MODEL}_P=${MOMENTUM}" --png
        $PLOT "TOTAL_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat"     "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"     "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"     "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"     --title "TOTAL ENERGY: ${MODEL}"         --xlabel "Time (a.u.)" --ylabel "Energy (a.u.)"     --legend "EXACT" "FSSH" "LZSH" "MASH" --output "TOTAL_ENERGY_${MODEL}_P=${MOMENTUM}"     --png
        $PLOT "POSITION_${MODEL}_P=${MOMENTUM}_EXACT.mat"         "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"         "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"         "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"         --title "POSITION: ${MODEL}"             --xlabel "Time (a.u.)" --ylabel "Coordinate (a.u.)" --legend "EXACT" "FSSH" "LZSH" "MASH" --output "POSITION_${MODEL}_P=${MOMENTUM}"         --png
        $PLOT "MOMENTUM_${MODEL}_P=${MOMENTUM}_EXACT.mat"         "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"         "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"         "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"         --title "MOMENTUM: ${MODEL}"             --xlabel "Time (a.u.)" --ylabel "Momentum (a.u.)"   --legend "EXACT" "FSSH" "LZSH" "MASH" --output "MOMENTUM_${MODEL}_P=${MOMENTUM}"         --png

        $PLOT "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"             "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"             --title "ADIABATIC POPULATION: ${MODEL}"     --xlabel "Time (a.u.)" --ylabel "State Population" --legend "${LEGEND_FSSH[@]}"      --output "POPULATION_FSSH_BAECKAN_${MODEL}_P=${MOMENTUM}" --png
        $PLOT "TDC_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat:${INDICES_TDC::-1}" "TDC_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat:${INDICES_TDC::-1}" --title "TIME DERIVATIVE COUPLING: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "TDC (??)"         --legend "${LEGEND_TDC[@]}"       --output "TDC_FSSH_BAECKAN_${MODEL}_P=${MOMENTUM}"        --png
        $PLOT "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"             "COEFFICIENT_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"            --title "ADIABATIC POPULATION: ${MODEL}"     --xlabel "Time (a.u.)" --ylabel "State Population" --legend "${LEGEND_FSSH_COEF[@]}" --output "COEFFICIENT_FSSH_${MODEL}_P=${MOMENTUM}"        --png
        $PLOT "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"             "COEFFICIENT_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"            --title "ADIABATIC POPULATION: ${MODEL}"     --xlabel "Time (a.u.)" --ylabel "State Population" --legend "${LEGEND_KTSH_COEF[@]}" --output "COEFFICIENT_KTSH_${MODEL}_P=${MOMENTUM}"        --png

        montage "POTENTIAL_ADIABATIC_${MODEL}.png" "POPULATION_${MODEL}_P=${MOMENTUM}.png"              "POSITION_${MODEL}_P=${MOMENTUM}.png"         "MOMENTUM_${MODEL}_P=${MOMENTUM}.png"                                                        -mode concatenate -tile x1 "TRAJECTORY_GENERAL_${MODEL}_P=${MOMENTUM}.png"
        montage "POTENTIAL_ADIABATIC_${MODEL}.png" "KINETIC_ENERGY_${MODEL}_P=${MOMENTUM}.png"          "POTENTIAL_ENERGY_${MODEL}_P=${MOMENTUM}.png" "TOTAL_ENERGY_${MODEL}_P=${MOMENTUM}.png"                                                    -mode concatenate -tile x1 "TRAJECTORY_ENERGY_${MODEL}_P=${MOMENTUM}.png"
        montage "POTENTIAL_ADIABATIC_${MODEL}.png" "POPULATION_FSSH_BAECKAN_${MODEL}_P=${MOMENTUM}.png" "TDC_FSSH_BAECKAN_${MODEL}_P=${MOMENTUM}.png" "COEFFICIENT_FSSH_${MODEL}_P=${MOMENTUM}.png" "COEFFICIENT_KTSH_${MODEL}_P=${MOMENTUM}.png" -mode concatenate -tile x1 "TRAJECTORY_COEFFICIENT_${MODEL}_P=${MOMENTUM}.png"
    done

    $PLOT "FINAL_POPULATION_P_${MODEL}.mat" --title "MOMENTUM DEPENDENT IS POPULATION: ${MODEL}" --xlabel "Momentum (a.u.)" --ylabel "Initial State Population" --legend "EXACT" "FSSH" "KTSH" "LZSH" "MASH" --output "FINAL_POPULATION_P_${MODEL}" --png

    montage "POTENTIAL_DIABATIC_${MODEL}.png"  "POTENTIAL_ADIABATIC_${MODEL}.png" -mode concatenate -tile x1        "POTENTIAL_${MODEL}.png"
    montage "POTENTIAL_ADIABATIC_${MODEL}.png" "FINAL_POPULATION_P_${MODEL}.png"  -mode concatenate -tile x1 "FINAL_POPULATION_${MODEL}.png"
done

FILES=(${MODELS[@]/#/POTENTIAL_});        montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x        "POTENTIAL.png"
FILES=(${MODELS[@]/#/FINAL_POPULATION_}); montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x "FINAL_POPULATION.png"

for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORY_GENERAL_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "TRAJECTORY_GENERAL_P=${MOMENTUM}.png"
done

for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORY_ENERGY_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "TRAJECTORY_ENERGY_P=${MOMENTUM}.png"
done

for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORY_COEFFICIENT_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "TRAJECTORY_COEFFICIENT_P=${MOMENTUM}.png"
done

for MODEL in ${MODELS[@]}; do rm *${MODEL}*.png; done
