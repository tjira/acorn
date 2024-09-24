#!/bin/bash

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
        "energy"              : true
    }
}
EOM

read -r -d '' TEMPLATE_SH_DYN <<- EOM
{
    "adiabatic" : true,
    "iterations" : 5000,
    "time_step" : 1,
    "trajectories" : 1000,
    "log_interval" : 1000,
    "data_export" : {
        "diabatic_population" : true, "adiabatic_population" : true,
        "position"            : true, "momentum"             : true,
        "energy"              : true
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

# define the model names
MODELS=("TULLY_1" "TULLY_2" "DS_1" "DS_2" "TS_1" "TS_2" "TS_3" "TS_4")

# generate momenta
MOMENTA=($(seq 10.0 1.0 50.0))

# loop over the models
for MODEL in ${MODELS[@]}; do

    # create the last population file
    echo "${#MOMENTA[@]} 4" > "${MODEL}_FINAL_POPULATIONS.mat"

    # loop over the momentum
    for MOMENTUM in ${MOMENTA[@]}; do

        # create the json files
        echo "{}" > exact_${MODEL,,}_P=${MOMENTUM}.json && echo "{}" > fssh_${MODEL,,}_P=${MOMENTUM}.json && echo "{}" > lzsh_${MODEL,,}_P=${MOMENTUM}.json

        # fill the json files with the templated wavefunctions
        jq --arg wavefunction "${TEMPLATE_WAVEFUNCTION}" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg wavefunction "${TEMPLATE_WAVEFUNCTION}" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg wavefunction "${TEMPLATE_WAVEFUNCTION}" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # fill the json files with the templated dynamics
        jq --arg dynamics "${TEMPLATE_EXACT_DYN}" '. |= . + {"quantum_dynamics"   : ($dynamics | fromjson)}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg dynamics "${TEMPLATE_SH_DYN}"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg dynamics "${TEMPLATE_SH_DYN}"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # fill the json files with the potential
        jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.quantum_dynamics   |= . + {"potential" : ($potential | fromjson)}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.classical_dynamics |= . + {"potential" : ($potential | fromjson)}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
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
        jq --arg guess "${GUESS}" '.wavefunction |= . + {"guess" : ($guess | fromjson)}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # set the momentum
        jq --arg momentum "${MOMENTUM}" '.wavefunction |= . + {"momentum" : ($momentum | fromjson)}' "exact_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg momentum "${MOMENTUM}" '.wavefunction |= . + {"momentum" : ($momentum | fromjson)}' "fssh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq --arg momentum "${MOMENTUM}" '.wavefunction |= . + {"momentum" : ($momentum | fromjson)}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json"  > temp.json && mv temp.json  "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # fill the json files with the surface hopping type
        jq '.classical_dynamics |= . + {"surface_hopping" : {"type" : "fewest-switches"}}' "fssh_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "fssh_${MODEL,,}_P=${MOMENTUM}.json"
        jq '.classical_dynamics |= . + {"surface_hopping" : {"type" : "landau-zener"   }}' "lzsh_${MODEL,,}_P=${MOMENTUM}.json" > temp.json && mv temp.json "lzsh_${MODEL,,}_P=${MOMENTUM}.json"

        # run the dynamics
        acorn -i "exact_${MODEL,,}_P=${MOMENTUM}.json" "fssh_${MODEL,,}_P=${MOMENTUM}.json" "lzsh_${MODEL,,}_P=${MOMENTUM}.json" -n 1

        # get the populations at the last time step
        POP_EXACT=$(tail -n 1 POPULATION_ADIABATIC_EXACT_REAL_1.mat | awk -v i=$IS '{print $(i+1)}');
        POP_FSSH=$( tail -n 1 POPULATION_ADIABATIC_FS-ADIABATIC.mat | awk -v i=$IS '{print $(i+1)}');
        POP_LZSH=$( tail -n 1 POPULATION_ADIABATIC_LZ-ADIABATIC.mat | awk -v i=$IS '{print $(i+1)}')

        # append the populations to the file
        echo "${MOMENTUM} ${POP_EXACT} ${POP_FSSH} ${POP_LZSH}" >> "${MODEL}_FINAL_POPULATIONS.mat"

        # fill the legend array for the population
        for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} EXACT"); done
        for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} FSSH" ); done
        for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} LZSH" ); done

        # fill the legend array for the potentials
        for (( i=0; i<$STATES; i++ )); do LEGEND_POT_DIA+=("S${i}"); LEGEND_POT_ADIA+=("S${i}"); done

        # get the indices on the diabatic potential diagonal
        DIA_INDICES=""; ADIA_INDICES=""; for (( i=0; i<$STATES; i++ )); do DIA_INDICES+=$(echo "$i*$STATES+$i" | bc -l)","; ADIA_INDICES+="$i,"; done

        # plot the population
        plot-1d.py POPULATION_ADIABATIC_EXACT_REAL_1.mat POPULATION_ADIABATIC_FS-ADIABATIC.mat POPULATION_ADIABATIC_LZ-ADIABATIC.mat --legend "${LEGEND_POP[@]}" --title "ADIABATIC POPULATION: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "State Population" --output "POPULATION_ADIABATIC_${MODEL}_P=${MOMENTUM}" --png

        # plot the potentials
        plot-1d.py "POTENTIAL_ADIABATIC.mat:${ADIA_INDICES::-1}" --legend "${LEGEND_POT_ADIA[@]}" --title "ADIABATIC POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)" --output "POTENTIAL_ADIABATIC_${MODEL}" --png
        plot-1d.py "POTENTIAL_DIABATIC.mat:${DIA_INDICES::-1}"   --legend "${LEGEND_POT_DIA[@]}"  --title "DIABATIC_POTENTIAL: ${MODEL}" --xlabel "Coordinate (a.u.)" --ylabel "Energy (a.u.)"  --output "POTENTIAL_DIABATIC_${MODEL}"  --png

        # average the position and momentum for classically propagated trajectories
        awk '{printf((NR == 1) ? "%d " : "%.2f ", $1); sum=0; for(i=2;i<=NF;i++) sum+=$i; print (NR == 1) ? 1 : sum/NF}' POSITION_LZ-ADIABATIC.mat > POSITION_LZ.mat
        awk '{printf((NR == 1) ? "%d " : "%.2f ", $1); sum=0; for(i=2;i<=NF;i++) sum+=$i; print (NR == 1) ? 1 : sum/NF}' MOMENTUM_LZ-ADIABATIC.mat > MOMENTUM_LZ.mat
        awk '{printf((NR == 1) ? "%d " : "%.2f ", $1); sum=0; for(i=2;i<=NF;i++) sum+=$i; print (NR == 1) ? 1 : sum/NF}' ENERGY_LZ-ADIABATIC.mat   >   ENERGY_LZ.mat
        awk '{printf((NR == 1) ? "%d " : "%.2f ", $1); sum=0; for(i=2;i<=NF;i++) sum+=$i; print (NR == 1) ? 1 : sum/NF}' POSITION_FS-ADIABATIC.mat > POSITION_FS.mat
        awk '{printf((NR == 1) ? "%d " : "%.2f ", $1); sum=0; for(i=2;i<=NF;i++) sum+=$i; print (NR == 1) ? 1 : sum/NF}' MOMENTUM_FS-ADIABATIC.mat > MOMENTUM_FS.mat
        awk '{printf((NR == 1) ? "%d " : "%.2f ", $1); sum=0; for(i=2;i<=NF;i++) sum+=$i; print (NR == 1) ? 1 : sum/NF}' ENERGY_FS-ADIABATIC.mat   >   ENERGY_FS.mat

        # plot the position and momentum
        plot-1d.py POSITION_EXACT_REAL_1.mat POSITION_FS.mat POSITION_LZ.mat --legend "EXACT" "FSSH" "LZSH" --title "POSITION: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "Position (a.u.)" --output "POSITION_${MODEL}_P=${MOMENTUM}" --png
        plot-1d.py MOMENTUM_EXACT_REAL_1.mat MOMENTUM_FS.mat MOMENTUM_LZ.mat --legend "EXACT" "FSSH" "LZSH" --title "MOMENTUM: ${MODEL}" --xlabel "Time (a.u.)" --ylabel "Momentum (a.u.)" --output "MOMENTUM_${MODEL}_P=${MOMENTUM}" --png
        plot-1d.py ENERGY_EXACT_REAL_1.mat   ENERGY_FS.mat   ENERGY_LZ.mat   --legend "EXACT" "FSSH" "LZSH" --title "ENERGY: ${MODEL}"   --xlabel "Time (a.u.)" --ylabel "ENERGY (a.u.)"   --output "ENERGY_${MODEL}_P=${MOMENTUM}"   --png

        # make the trajectory analysis image
        montage "POTENTIAL_ADIABATIC_${MODEL}.png" "POPULATION_ADIABATIC_${MODEL}_P=${MOMENTUM}.png" "POSITION_${MODEL}_P=${MOMENTUM}.png" "MOMENTUM_${MODEL}_P=${MOMENTUM}.png" "ENERGY_${MODEL}_P=${MOMENTUM}.png" -mode concatenate -tile x1 "TRAJECTORIES_${MODEL}_P=${MOMENTUM}".png
    done

    # plot the population dependence on momentum
    plot-1d.py "${MODEL}_FINAL_POPULATIONS.mat" --legend "EXACT" "FSSH" "LZSH" --title "FINAL POPULATION: ${MODEL}" --xlabel "Initial Momentum (a.u.)" --ylabel "Final Population of the Initial State (S$IS)" --output "${MODEL}_FINAL_POPULATIONS" --png

    # create the potential plot
    montage "POTENTIAL_DIABATIC_${MODEL}.png" "POTENTIAL_ADIABATIC_${MODEL}.png" -mode concatenate -tile x1 "POTENTIALS_${MODEL}.png"

    # make the comparison image
    montage "POTENTIAL_ADIABATIC_${MODEL}.png" "${MODEL}_FINAL_POPULATIONS.png" -mode concatenate -tile x1 "COMPARISONS_${MODEL}.png"
done

# montage all the images together
for MOMENTUM in ${MOMENTA[@]}; do
    FILES=(${MODELS[@]/#/TRAJECTORIES_}); montage "${FILES[@]/%/_P=${MOMENTUM}.png}" -mode concatenate -tile 1x "TRAJECTORIES_P=${MOMENTUM}.png"
done

# montage all the potential images together
FILES=(${MODELS[@]/#/POTENTIALS_}); montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x POTENTIALS.png

# montage all the comparison images together
FILES=(${MODELS[@]/#/COMPARISONS_}); montage "${FILES[@]/%/.png}" -mode concatenate -tile 1x COMPARISONS.png

# remove the intermediate files
mkdir -p .temp && mv COMPARISONS.png POTENTIALS.png TRAJECTORIES_P* .temp && rm *.json *.mat *.png && mv .temp/* . && rm -r .temp
