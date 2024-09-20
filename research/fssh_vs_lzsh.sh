#!/bin/bash

# ======================================================================================================================================================================================================
# START TEMPLATES
# ======================================================================================================================================================================================================

read -r -d '' TEMPLATE_WAVEFUNCTION <<- EOM
{
    "dimension" : 1,
    "grid_limits" : [-24, 40],
    "grid_points" : 1024,
    "mass" : 2000,
    "momentum" : 15
}
EOM

read -r -d '' TEMPLATE_EXACT_DYN <<- EOM
{
    "real" : 1,
    "iterations" : 350,
    "time_step" : 10,
    "data_export" : {
        "diabatic_population" : true, "adiabatic_population" : true,
        "diabatic_potential"  : true, "adiabatic_potential"  : true
    }
}
EOM

read -r -d '' TEMPLATE_SH_DYN <<- EOM
{
    "adiabatic" : true,
    "iterations" : 3500,
    "time_step" : 1,
    "trajectories" : 1000,
    "log_interval" : 500,
    "data_export" : {
        "diabatic_population" : true, "adiabatic_population" : true
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

MODELS=("TULLY_1" "TULLY_2" "DS_1" "DS_2" "TS_1" "TS_2" "TS_3" "TS_4")

# loop over the models
for MODEL in ${MODELS[@]}; do

    # create the json files
    echo "{}" > exact_${MODEL,,}.json && echo "{}" > fssh_${MODEL,,}.json && echo "{}" > lzsh_${MODEL,,}.json

    # fill the json files with the templated wavefunctions
    jq --arg wavefunction "$TEMPLATE_WAVEFUNCTION" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' exact_${MODEL,,}.json > temp.json && mv temp.json exact_${MODEL,,}.json
    jq --arg wavefunction "$TEMPLATE_WAVEFUNCTION" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' fssh_${MODEL,,}.json  > temp.json && mv temp.json  fssh_${MODEL,,}.json
    jq --arg wavefunction "$TEMPLATE_WAVEFUNCTION" '. |= . + {"wavefunction" : ($wavefunction | fromjson)}' lzsh_${MODEL,,}.json  > temp.json && mv temp.json  lzsh_${MODEL,,}.json

    # fill the json files with the templated dynamics
    jq --arg dynamics "$TEMPLATE_EXACT_DYN" '. |= . + {"quantum_dynamics"   : ($dynamics | fromjson)}' exact_${MODEL,,}.json > temp.json && mv temp.json exact_${MODEL,,}.json
    jq --arg dynamics "$TEMPLATE_SH_DYN"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}' fssh_${MODEL,,}.json  > temp.json && mv temp.json  fssh_${MODEL,,}.json
    jq --arg dynamics "$TEMPLATE_SH_DYN"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}' lzsh_${MODEL,,}.json  > temp.json && mv temp.json  lzsh_${MODEL,,}.json

    # fill the json files with the potential
    jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.quantum_dynamics   |= . + {"potential" : ($potential | fromjson)}' exact_${MODEL,,}.json > temp.json && mv temp.json exact_${MODEL,,}.json
    jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.classical_dynamics |= . + {"potential" : ($potential | fromjson)}' fssh_${MODEL,,}.json  > temp.json && mv temp.json  fssh_${MODEL,,}.json
    jq --arg potential "$(eval echo \$"POTENTIAL_${MODEL}")" '.classical_dynamics |= . + {"potential" : ($potential | fromjson)}' lzsh_${MODEL,,}.json  > temp.json && mv temp.json  lzsh_${MODEL,,}.json

    # get the number of states and create the guess and legend array
    STATES=$(jq '.quantum_dynamics.potential | length' exact_${MODEL,,}.json); GUESS='['; LEGEND_POP=(); LEGEND_POT_DIA=(); LEGEND_POT_ADIA=()

    # fil the initial guess for the wavefunction
    for (( i=0; i<$(eval echo \$"IS_${MODEL}"); i++ )); do GUESS+='"0","0",'; done; GUESS+='"exp(-(x+15)^2)","0",'; for (( i=0; i<$STATES-$(eval echo \$"IS_${MODEL}")-1; i++ )); do GUESS+='"0","0",'; done;GUESS=${GUESS::-1}]

    # add the guess to the json files
    jq --arg guess "$GUESS" '.wavefunction |= . + {"guess" : ($guess | fromjson)}' exact_${MODEL,,}.json > temp.json && mv temp.json exact_${MODEL,,}.json
    jq --arg guess "$GUESS" '.wavefunction |= . + {"guess" : ($guess | fromjson)}' fssh_${MODEL,,}.json  > temp.json && mv temp.json  fssh_${MODEL,,}.json
    jq --arg guess "$GUESS" '.wavefunction |= . + {"guess" : ($guess | fromjson)}' lzsh_${MODEL,,}.json  > temp.json && mv temp.json  lzsh_${MODEL,,}.json

    # fill the json files with the surface hopping type
    jq '.classical_dynamics |= . + {"surface_hopping" : {"type" : "fewest-switches"}}' fssh_${MODEL,,}.json > temp.json && mv temp.json fssh_${MODEL,,}.json
    jq '.classical_dynamics |= . + {"surface_hopping" : {"type" : "landau-zener"   }}' lzsh_${MODEL,,}.json > temp.json && mv temp.json lzsh_${MODEL,,}.json

    # run the dynamics
    acorn -i exact_${MODEL,,}.json fssh_${MODEL,,}.json lzsh_${MODEL,,}.json

    # fill the legend array for the population
    for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} EXACT"); done
    for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} FSSH" ); done
    for (( i=0; i<$STATES; i++ )); do LEGEND_POP+=("S${i} LZSH" ); done

    # fill the legend array for the potentials
    for (( i=0; i<$STATES; i++ )); do LEGEND_POT_DIA+=("S${i}"); LEGEND_POT_ADIA+=("S${i}"); done

    # get the indices on the diabatic potential diagonal
    DIA_INDICES=""; ADIA_INDICES=""; for (( i=0; i<$STATES; i++ )); do DIA_INDICES+=$(echo "$i*$STATES+$i" | bc -l)","; ADIA_INDICES+="$i,"; done

    # plot the population
    plot-1d.py POPULATION_ADIABATIC_EXACT_REAL_1.mat POPULATION_ADIABATIC_FS-ADIABATIC.mat POPULATION_ADIABATIC_LZ-ADIABATIC.mat --legend "${LEGEND_POP[@]}" --title "ADIABATIC POPULATION: $MODEL" --output POPULATION_ADIABATIC_$MODEL --png

    # plot the potentials
    plot-1d.py POTENTIAL_ADIABATIC.mat":"${ADIA_INDICES::-1} --legend "${LEGEND_POT_ADIA[@]}" --title "ADIABATIC POTENTIAL: $MODEL" --output POTENTIAL_ADIABATIC_$MODEL --png
    plot-1d.py POTENTIAL_DIABATIC.mat":"${DIA_INDICES::-1}   --legend "${LEGEND_POT_DIA[@]}"  --title "DIABATIC_POTENTIAL: $MODEL"  --output POTENTIAL_DIABATIC_$MODEL  --png

    # join the potentials and the population in a row
    montage POTENTIAL_DIABATIC_$MODEL.png POTENTIAL_ADIABATIC_$MODEL.png POPULATION_ADIABATIC_$MODEL.png -mode concatenate -tile x1 $MODEL.png && rm *.json *.mat
done

# define png files to remove
PNG_REMOVE=(${MODELS[@]/#/\*}); PNG_REMOVE=(${PNG_REMOVE[@]/%/.png})

# create the final plot and remove the intermediate files
montage "${MODELS[@]/%/.png}" -mode concatenate -tile 1x dynamics.png && rm "${PNG_REMOVE[@]}"
