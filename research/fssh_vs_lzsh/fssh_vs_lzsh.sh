#!/bin/bash

# ======================================================================================================================================================================================================
# VARIABLES
# ======================================================================================================================================================================================================

CORES=1; PSTART=10.0; PSTEP=10.0; PEND=50.0; TRAJS=1000; LOG_INTERVAL_STEP=1000; LOG_INTERVAL_TRAJ=100; CLEAN=1

MODELS=("TULLY_1" "TULLY_2" "DS_1" "DS_2" "TS_1" "TS_2" "TS_3" "TS_4")

# ======================================================================================================================================================================================================
# EXECUTABLE PATHS
# ======================================================================================================================================================================================================

ACORN=./../../bin/acorn; PLOT_1D=./../../script/plot/plot-1d.py

# ======================================================================================================================================================================================================
# START TEMPLATES
# ======================================================================================================================================================================================================

read -r -d '' TEMPLATE_WAVEFUNCTION <<- EOM
{
    "dimension" : 1,
    "grid_limits" : [-48, 144],
    "grid_points" : 4096,
    "mass" : 2000
}
EOM

read -r -d '' TEMPLATE_EXACT_DYN <<- EOM
{
    "real" : 1,
    "iterations" : 500,
    "time_step" : 10,
    "data_export" : {
        "diabatic_population" : true, "adiabatic_population" : true,
        "diabatic_potential"  : true, "adiabatic_potential"  : true,
        "position"            : true, "momentum"             : true,
        "total_energy"        : true, "potential_energy"     : true,
        "kinetic_energy"      : true
    }
}
EOM

read -r -d '' TEMPLATE_SH_DYN <<- EOM
{
    "adiabatic" : true,
    "iterations" : 5000,
    "time_step" : 1,
    "trajectories" : $TRAJS,
    "log_interval_step" : $LOG_INTERVAL_STEP,
    "log_interval_traj" : $LOG_INTERVAL_TRAJ,
    "data_export" : {
        "diabatic_population" : true, "adiabatic_population" : true,
        "position_mean"       : true, "momentum_mean"        : true,
        "total_energy_mean"   : true, "potential_energy_mean": true,
        "kinetic_energy_mean" : true, "hopping_geometry"     : true,
        "hopping_time"        : true
    }
}
EOM

# ======================================================================================================================================================================================================
# END TEMPLATES
# ======================================================================================================================================================================================================
# START POTENTIALS
# ======================================================================================================================================================================================================

read -r -d '' POTENTIAL_TULLY_1 <<- EOM
[
    ["if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))", "0.005*exp(-x^2)"                                   ],
    ["0.005*exp(-x^2)",                                   "-if(x>0,0.01*(1-exp(-1.6*x)),-0.01*(1-exp(1.6*x)))"]
]
EOM
IS_TULLY_1=1

read -r -d '' POTENTIAL_TULLY_2 <<- EOM
[
    ["0",                    "0.015*exp(-0.06*x^2)"    ],
    ["0.015*exp(-0.06*x^2)", "-0.1*exp(-0.28*x^2)+0.05"]
]
EOM
IS_TULLY_2=1

read -r -d '' POTENTIAL_DS_1 <<- EOM
[
    ["0.001*x",               "0.001*exp(-0.05*x^2)"],
    ["0.001*exp(-0.05*x^2)",  "-0.001*x"            ]
]
EOM
IS_DS_1=1

read -r -d '' POTENTIAL_DS_2 <<- EOM
[
    ["0.01*tanh(0.6*x)", "0.001*exp(-x^2)"  ],
    ["0.001*exp(-x^2)",  "-0.01*tanh(0.6*x)"]
]
EOM
IS_DS_2=1

read -r -d '' POTENTIAL_TS_1 <<- EOM
[
    ["0.001*(x-5)",              "0.001*exp(-0.01*(x-5)^2)", "0.002*exp(-0.01*x^2)"    ],
    ["0.001*exp(-0.01*(x-5)^2)", "0",                        "0.001*exp(-0.01*(x+5)^2)"],
    ["0.002*exp(-0.01*x^2)",     "0.001*exp(-0.01*(x+5)^2)", "-0.001*(x+5)"            ]
]
EOM
IS_TS_1=2

read -r -d '' POTENTIAL_TS_2 <<- EOM
[
    ["0.03*(tanh(1.6*x)+tanh(1.6*(x+7)))", "0.005*exp(-x^2)",                     "0.005*exp(-(x+7)^2)"                    ],
    ["0.005*exp(-x^2)",                    "-0.03*(tanh(1.6*x)+tanh(1.6*(x-7)))", "0.005*exp(-(x-7)^2)"                    ],
    ["0.005*exp(-(x+7)^2)",                "0.005*exp(-(x-7)^2)",                 "-0.03*(tanh(1.6*(x+7))-tanh(1.6*(x-7)))"]
]
EOM
IS_TS_2=1

read -r -d '' POTENTIAL_TS_3 <<- EOM
[
    ["0.001*x",              "0.001*exp(-0.01*x^2)", "0.002*exp(-0.01*x^2)"],
    ["0.001*exp(-0.01*x^2)", "0",                    "0.001*exp(-0.01*x^2)"],
    ["0.002*exp(-0.01*x^2)", "0.001*exp(-0.01*x^2)", "-0.001*x"            ]
]
EOM
IS_TS_3=2

read -r -d '' POTENTIAL_TS_4 <<- EOM
[
    ["0.01*tanh(0.5*x)", "0.001*exp(-x^2)", "0.002*exp(-x^2)"  ],
    ["0.001*exp(-x^2)",  "0",               "0.001*exp(-x^2)"  ],
    ["0.002*exp(-x^2)",  "0.001*exp(-x^2)", "-0.01*tanh(0.5*x)"]
]
EOM
IS_TS_4=2

# ======================================================================================================================================================================================================
# END POTENTIALS
# ======================================================================================================================================================================================================

# remove all the previous files
[[ $CLEAN -eq 1 ]] && rm -f *.json *.mat *.png

# generate momenta
(( $(echo "$PSTEP <= 0" | bc -l) )) && MOMENTA=($PSTART) || MOMENTA=($(seq $PSTART $PSTEP $PEND))

# loop over the models
for MODEL in ${MODELS[@]}; do

    # create the last population file
    echo "${#MOMENTA[@]} 4" > "${MODEL}_FINAL_POPULATIONS.mat"

    # loop over the momentum
    for MOMENTUM in ${MOMENTA[@]}; do

        # create the json files
        echo "{}" > exact_${MODEL,,}_P=${MOMENTUM}.json && echo "{}" > fssh_${MODEL,,}_P=${MOMENTUM}.json && echo "{}" > kfssh_${MODEL,,}_P=${MOMENTUM}.json && echo "{}" > lzsh_${MODEL,,}_P=${MOMENTUM}.json

        # fill the json files with the templated wavefunctions
        jq --arg wavefunction "${TEMPLATE_WAVEFUNCTION}" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg wavefunction "${TEMPLATE_WAVEFUNCTION}" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg wavefunction "${TEMPLATE_WAVEFUNCTION}" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' "kfssh_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg wavefunction "${TEMPLATE_WAVEFUNCTION}" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # fill the json files with the templated dynamics
        jq --arg dynamics "${TEMPLATE_EXACT_DYN}" '. |= . + {"quantum_dynamics"   : ($dynamics | fromjson)}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg dynamics "${TEMPLATE_SH_DYN}"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg dynamics "${TEMPLATE_SH_DYN}"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}' "kfssh_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg dynamics "${TEMPLATE_SH_DYN}"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # fill the json files with the potential
        jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.quantum_dynamics   |= . + {"potential" : ($potential | fromjson)}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.classical_dynamics |= . + {"potential" : ($potential | fromjson)}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.classical_dynamics |= . + {"potential" : ($potential | fromjson)}' "kfssh_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.classical_dynamics |= . + {"potential" : ($potential | fromjson)}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # get the number of states and create the guess and legend array
        STATES=$(jq '.quantum_dynamics.potential | length' "exact_${MODEL,,}_P=${MOMENTUM}.json"); GUESS='['; LEGEND_POP=(); LEGEND_POT_DIA=(); LEGEND_POT_ADIA=()

        # extract the initial state
        IS=$(eval echo \$"IS_${MODEL}")

        # fil the initial guess for the wavefunction
        for (( i=0; i<$IS; i++ )); do GUESS+='"0","0",'; done; GUESS+='"exp(-(x+15)^2)","0",'; for (( i=0; i<$STATES-$IS-1; i++ )); do GUESS+='"0","0",'; done;GUESS=${GUESS::-1}]

        # add the guess to the json files
        jq --arg guess "${GUESS}" '.wavefunction |= . + {"guess" : ($guess | fromjson)}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg guess "${GUESS}" '.wavefunction |= . + {"guess" : ($guess | fromjson)}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg guess "${GUESS}" '.wavefunction |= . + {"guess" : ($guess | fromjson)}' "kfssh_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg guess "${GUESS}" '.wavefunction |= . + {"guess" : ($guess | fromjson)}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # set the momentum
        jq --arg momentum "${MOMENTUM}" '.wavefunction |= . + {"momentum" : [($momentum | fromjson)]}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg momentum "${MOMENTUM}" '.wavefunction |= . + {"momentum" : [($momentum | fromjson)]}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg momentum "${MOMENTUM}" '.wavefunction |= . + {"momentum" : [($momentum | fromjson)]}' "kfssh_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg momentum "${MOMENTUM}" '.wavefunction |= . + {"momentum" : [($momentum | fromjson)]}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # fill the json files with the surface hopping type
        jq '.classical_dynamics |= . + {"surface_hopping" : {"type" : "fewest-switches", "kappa" : false}}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq '.classical_dynamics |= . + {"surface_hopping" : {"type" : "fewest-switches", "kappa" : true }}' "kfssh_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq '.classical_dynamics |= . + {"surface_hopping" : {"type" : "landau-zener"                    }}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # run the dynamics and delete the inputs
        $ACORN -i "exact_${MODEL,,}_P=${MOMENTUM}.json" "fssh_${MODEL,,}_P=${MOMENTUM}.json" "kfssh_${MODEL,,}_P=${MOMENTUM}.json" "lzsh_${MODEL,,}_P=${MOMENTUM}.json" -n $CORES && rm *.json

        # get the populations at the last time step
        POP_EXACT=$(tail -n 1 POPULATION-ADIABATIC_EXACT-REAL_1.mat  | awk -v i=$IS '{print $(i+1)}');
        POP_FSSH=$( tail -n 1 POPULATION-ADIABATIC_FS-ADIABATIC.mat  | awk -v i=$IS '{print $(i+1)}');
        POP_KFSSH=$(tail -n 1 POPULATION-ADIABATIC_KFS-ADIABATIC.mat | awk -v i=$IS '{print $(i+1)}');
        POP_LZSH=$( tail -n 1 POPULATION-ADIABATIC_LZ-ADIABATIC.mat  | awk -v i=$IS '{print $(i+1)}')

        # append the populations to the file
        echo "${MOMENTUM} ${POP_EXACT} ${POP_FSSH} ${POP_KFSSH} ${POP_LZSH}" >> "${MODEL}_FINAL_POPULATIONS.mat"

        # fill the legend array for the population
        for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} EXACT"); done
        for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} FSSH" ); done
        for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} KFSSH"); done
        for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} LZSH" ); done

        # fill the legend array for the potentials
        for (( i=0; i<$STATES; i++ )); do LEGEND_POT_DIA+=("S${i}"); LEGEND_POT_ADIA+=("S${i}"); done

        # get the indices on the diabatic potential diagonal
        DIA_INDICES=""; ADIA_INDICES=""; for (( i=0; i<$STATES; i++ )); do DIA_INDICES+=$(echo "$i*$STATES+$i" | bc -l)","; ADIA_INDICES+="$i,"; done

        # plot the population
        $PLOT_1D POPULATION-ADIABATIC_EXACT-REAL_1.mat POPULATION-ADIABATIC_FS-ADIABATIC.mat POPULATION-ADIABATIC_KFS-ADIABATIC.mat POPULATION-ADIABATIC_LZ-ADIABATIC.mat --legend "${LEGEND_POP[@]}" --title "ADIABATIC POPULATION: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "State Population" --output "POPULATION-ADIABATIC_${MODEL}_P=${MOMENTUM}" --png

        # plot the potentials
        $PLOT_1D "POTENTIAL-ADIABATIC.mat:${ADIA_INDICES::-1}" --legend "${LEGEND_POT_ADIA[@]}" --title "ADIABATIC POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --output "POTENTIAL-ADIABATIC_${MODEL}" --png
        $PLOT_1D "POTENTIAL-DIABATIC.mat:${DIA_INDICES::-1}"   --legend "${LEGEND_POT_DIA[@]}"  --title "DIABATIC POTENTIAL: ${MODEL}"  --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --output  "POTENTIAL-DIABATIC_${MODEL}" --png

        # plot the position and momentum
        $PLOT_1D POSITION_EXACT-REAL_1.mat         POSITION-MEAN_FS-ADIABATIC.mat         POSITION-MEAN_KFS-ADIABATIC.mat         POSITION-MEAN_LZ-ADIABATIC.mat         --legend "EXACT" "FSSH" "KFSSH" "LZSH" --title "POSITION: ${MODEL}"         --xlabel "Time (a.u.)" --ylabel "Position (a.u.)" --output "POSITION_${MODEL}_P=${MOMENTUM}"         --png
        $PLOT_1D MOMENTUM_EXACT-REAL_1.mat         MOMENTUM-MEAN_FS-ADIABATIC.mat         MOMENTUM-MEAN_KFS-ADIABATIC.mat         MOMENTUM-MEAN_LZ-ADIABATIC.mat         --legend "EXACT" "FSSH" "KFSSH" "LZSH" --title "MOMENTUM: ${MODEL}"         --xlabel "Time (a.u.)" --ylabel "Momentum (a.u.)" --output "MOMENTUM_${MODEL}_P=${MOMENTUM}"         --png
        $PLOT_1D TOTAL-ENERGY_EXACT-REAL_1.mat     TOTAL-ENERGY-MEAN_FS-ADIABATIC.mat     TOTAL-ENERGY-MEAN_KFS-ADIABATIC.mat     TOTAL-ENERGY-MEAN_LZ-ADIABATIC.mat     --legend "EXACT" "FSSH" "KFSSH" "LZSH" --title "TOTAL ENERGY: ${MODEL}"     --xlabel "Time (a.u.)" --ylabel "ENERGY (a.u.)"   --output "TOTAL-ENERGY_${MODEL}_P=${MOMENTUM}"     --png
        $PLOT_1D KINETIC-ENERGY_EXACT-REAL_1.mat   KINETIC-ENERGY-MEAN_FS-ADIABATIC.mat   KINETIC-ENERGY-MEAN_KFS-ADIABATIC.mat   KINETIC-ENERGY-MEAN_LZ-ADIABATIC.mat   --legend "EXACT" "FSSH" "KFSSH" "LZSH" --title "KINETIC ENERGY: ${MODEL}"   --xlabel "Time (a.u.)" --ylabel "ENERGY (a.u.)"   --output "KINETIC-ENERGY_${MODEL}_P=${MOMENTUM}"   --png
        $PLOT_1D POTENTIAL-ENERGY_EXACT-REAL_1.mat POTENTIAL-ENERGY-MEAN_FS-ADIABATIC.mat POTENTIAL-ENERGY-MEAN_KFS-ADIABATIC.mat POTENTIAL-ENERGY-MEAN_LZ-ADIABATIC.mat --legend "EXACT" "FSSH" "KFSSH" "LZSH" --title "POTENTIAL ENERGY: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "ENERGY (a.u.)"   --output "POTENTIAL-ENERGY_${MODEL}_P=${MOMENTUM}" --png

        # plot the hopping geometries
        $PLOT_1D HOPPING-GEOMETRIES_FS-ADIABATIC.mat HOPPING-GEOMETRIES_KFS-ADIABATIC.mat --bins 100 --legend "FSSH" "KFSSH" --title "HOPPING GEOMETRIES: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Relative Count" --output "HOPPING-GEOMETRIES_${MODEL}" --histogram --png
        $PLOT_1D HOPPING-TIMES_FS-ADIABATIC.mat           HOPPING-TIMES_KFS-ADIABATIC.mat --bins 100 --legend "FSSH" "KFSSH" --title "HOPPING TIMES: ${MODEL}"      --xlabel "Time (a.u.)"       --ylabel "Relative Count" --output "HOPPING-TIMES_${MODEL}"      --histogram --png

        # make the trajectory analysis image
        montage "POTENTIAL-ADIABATIC_${MODEL}.png" "POPULATION-ADIABATIC_${MODEL}_P=${MOMENTUM}.png" "POSITION_${MODEL}_P=${MOMENTUM}.png"         "MOMENTUM_${MODEL}_P=${MOMENTUM}.png"     -mode concatenate -tile x1 "TRAJECTORIES-GENERAL_${MODEL}_P=${MOMENTUM}.png"
        montage "POTENTIAL-ADIABATIC_${MODEL}.png" "KINETIC-ENERGY_${MODEL}_P=${MOMENTUM}.png"       "POTENTIAL-ENERGY_${MODEL}_P=${MOMENTUM}.png" "TOTAL-ENERGY_${MODEL}_P=${MOMENTUM}.png" -mode concatenate -tile x1 "TRAJECTORIES-ENERGETICS_${MODEL}_P=${MOMENTUM}.png"
        montage "POTENTIAL-ADIABATIC_${MODEL}.png" "HOPPING-GEOMETRIES_${MODEL}.png"                 "HOPPING-TIMES_${MODEL}.png"                                                            -mode concatenate -tile x1 "TRAJECTORIES-TRANSITIONS_${MODEL}_P=${MOMENTUM}.png"
    done

    # plot the population dependence on momentum
    $PLOT_1D "${MODEL}_FINAL_POPULATIONS.mat" --legend "EXACT" "FSSH" "KFSSH" "LZSH" --title "FINAL POPULATION: ${MODEL}" --xlabel "Initial Momentum (a.u.)" --ylabel "Final Population of the Initial State (S$IS)" --output "${MODEL}_FINAL_POPULATIONS" --png

    # make the comparison image
    montage "POTENTIAL-ADIABATIC_${MODEL}.png" "${MODEL}_FINAL_POPULATIONS.png" -mode concatenate -tile x1 "COMPARISONS_${MODEL}.png"

    # create the potential plot
    montage "POTENTIAL-DIABATIC_${MODEL}.png" "POTENTIAL-ADIABATIC_${MODEL}.png" -mode concatenate -tile x1 "POTENTIALS_${MODEL}.png"
    
    # remove every image except for comparisons, potentials and trajectories
    mkdir -p .temp && mv COMPARISONS_*.png POTENTIALS_*.png TRAJECTORIES-*.png .temp && rm *.mat *.png && mv .temp/* . && rm -r .temp
done

# montage all the images together and remove the individual ones
for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORIES-GENERAL_}   ); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "TRAJECTORIES-GENERAL_P=${MOMENTUM}.png"    && rm "${FILES[@]/%/_P=${MOMENTUM}.png}"
    FILES=(${MODELS[@]/#/TRAJECTORIES-ENERGETICS_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "TRAJECTORIES-ENERGETICS_P=${MOMENTUM}.png" && rm "${FILES[@]/%/_P=${MOMENTUM}.png}"
done

# montage all the potential images together and remove the individual ones
FILES=(${MODELS[@]/#/POTENTIALS_}); montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x POTENTIALS.png && rm "${FILES[@]/%/.png}"

# montage all the comparison images together and remove the individual ones
FILES=(${MODELS[@]/#/COMPARISONS_}); montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x COMPARISONS.png && rm "${FILES[@]/%/.png}"
