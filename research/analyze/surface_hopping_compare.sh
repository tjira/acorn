#!/bin/bash

MODELS=($(ls -l POTENTIAL_DIABATIC* | grep -oP "(?<=POTENTIAL_DIABATIC_)[^_]+_[^_]" | sort | uniq)); MOMENTA=($(ls -l POPULATION* | grep -oP "P=\K[0-9.]+" | sort | uniq))

for MODEL in ${MODELS[@]}; do

    NSTATE=$(tail -n 1 POPULATION_${MODEL}_* | tail -n 1 | awk '{print NF - 1}');

    LEGEND_POP=(); LEGEND_FSSH=(); LEGEND_FSSH_COEF=(); LEGEND_KTSH_COEF=(); LEGEND_POT=(); LEGEND_TDC=(); INDICES_POT=""; INDICES_TDC="";

    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} EXACT"); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} FSSH" ); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} KTSH" ); done
    for (( i=0; i<$NSTATE; i++ )); do LEGEND_POP+=("S${i} LZSH" ); done

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

    $PLOT  "POTENTIAL_DIABATIC_${MODEL}.mat:${INDICES_POT::-1}" --title  "DIABATIC POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --legend "${LEGEND_POT[@]}" --output  "POTENTIAL_DIABATIC_${MODEL}.png"
    $PLOT "POTENTIAL_ADIABATIC_${MODEL}.mat:${INDICES_POT::-1}" --title "ADIABATIC POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --legend "${LEGEND_POT[@]}" --output "POTENTIAL_ADIABATIC_${MODEL}.png"

    echo "${#MOMENTA[@]} 4" > "FINAL_POPULATION_P_${MODEL}.mat"

    for MOMENTUM in ${MOMENTA[@]}; do

        [[ -f "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      ]] && IS_COL_EXACT=$(head -n 2 "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  ]] && IS_COL_FSSH=$( head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"  ]] && IS_COL_KTSH=$( head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"  | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  ]] && IS_COL_LZSH=$( head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"  ]] && IS_COL_MASH=$( head -n 2 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"  | tail -n 1 | awk '{is=2; pop=$2; for(i=3; i<=NF; i++) if ($i > pop) {is=i; pop=$i;} print is - 2}')

        [[ -f "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      ]] && POP_EXACT=$(tail -n 1 "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      | awk -v i=$IS_COL_EXACT '{print $(i+2)}');
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  ]] && POP_FSSH=$( tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  | awk -v i=$IS_COL_FSSH  '{print $(i+2)}');
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"  ]] && POP_KTSH=$( tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"  | awk -v i=$IS_COL_KTSH  '{print $(i+2)}');
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  ]] && POP_LZSH=$( tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  | awk -v i=$IS_COL_LZSH  '{print $(i+2)}');
        [[ -f "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"  ]] && POP_MASH=$( tail -n 1 "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_MASH.mat"  | awk -v i=$IS_COL_LZSH  '{print $(i+2)}');

        echo "${MOMENTUM} ${POP_EXACT} ${POP_FSSH} ${POP_KTSH} ${POP_LZSH}" >> "FINAL_POPULATION_P_${MODEL}.mat"

        $PLOT "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat" "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat" "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat" "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat" --title "ADIABATIC POPULATION: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "State Population"  --legend "${LEGEND_POP[@]}" --output "POPULATION_${MODEL}_P=${MOMENTUM}.png"

        $PLOT "TDC_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat:${INDICES_TDC::-1}" "TDC_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat:${INDICES_TDC::-1}" --title "TIME DERIVATIVE COUPLING: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "TDC (??)"         --legend "${LEGEND_TDC[@]}"       --output "TDC_FSSH_BAECKAN_${MODEL}_P=${MOMENTUM}.png"
        $PLOT "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"             "COEFFICIENT_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"            --title "ADIABATIC POPULATION: ${MODEL}"     --xlabel "Time (a.u.)" --ylabel "State Population" --legend "${LEGEND_FSSH_COEF[@]}" --output "COEFFICIENT_FSSH_${MODEL}_P=${MOMENTUM}.png"
        $PLOT "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"             "COEFFICIENT_MEAN_${MODEL}_P=${MOMENTUM}_KTSH.mat"            --title "ADIABATIC POPULATION: ${MODEL}"     --xlabel "Time (a.u.)" --ylabel "State Population" --legend "${LEGEND_KTSH_COEF[@]}" --output "COEFFICIENT_KTSH_${MODEL}_P=${MOMENTUM}.png"

        montage "POTENTIAL_ADIABATIC_${MODEL}.png" "POPULATION_${MODEL}_P=${MOMENTUM}.png"       "TDC_FSSH_BAECKAN_${MODEL}_P=${MOMENTUM}.png" -mode concatenate -tile x1 "TRAJECTORY_POPULATIONS_${MODEL}_P=${MOMENTUM}.png"
        montage "POTENTIAL_ADIABATIC_${MODEL}.png" "COEFFICIENT_FSSH_${MODEL}_P=${MOMENTUM}.png" "COEFFICIENT_KTSH_${MODEL}_P=${MOMENTUM}.png" -mode concatenate -tile x1 "TRAJECTORY_COEFFICIENT_${MODEL}_P=${MOMENTUM}.png"
    done

    $PLOT "FINAL_POPULATION_P_${MODEL}.mat" --title "MOMENTUM DEPENDENT IS POPULATION: ${MODEL}" --xlabel "Momentum (a.u.)" --ylabel "Initial State Population" --legend "EXACT" "FSSH" "KTSH" "LZSH" "MASH" --output "FINAL_POPULATION_P_${MODEL}.png"

    montage "POTENTIAL_DIABATIC_${MODEL}.png"  "POTENTIAL_ADIABATIC_${MODEL}.png" -mode concatenate -tile x1 "POTENTIAL_${MODEL}.png"
    montage "POTENTIAL_ADIABATIC_${MODEL}.png" "FINAL_POPULATION_P_${MODEL}.png"  -mode concatenate -tile x1 "MOMDEPEND_${MODEL}.png"
done

FILES=(${MODELS[@]/#/POTENTIAL_}); montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x "POTENTIAL.png"
FILES=(${MODELS[@]/#/MOMDEPEND_}); montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x "MOMDEPEND.png"

for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORY_POPULATIONS_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "TRAJECTORY_POPULATIONS_P=${MOMENTUM}.png"
done

for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORY_COEFFICIENT_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "TRAJECTORY_COEFFICIENT_P=${MOMENTUM}.png"
done

for MODEL in ${MODELS[@]}; do rm *${MODEL}*.png; done
