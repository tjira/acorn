#!/bin/bash

# ======================================================================================================================================================================================================
# DEFAULT VARIABLES
# ======================================================================================================================================================================================================

CORES=1; PSTART=10.0; PSTEP=10.0; PEND=50.0; TRAJS=1000; CLEAN=1

MODELS=("tully1D_1" "tully1D_2" "tully1D_3" "doubleState1D_1" "doubleState1D_2" "tripleState1D_1" "tripleState1D_2" "tripleState1D_3")

# ======================================================================================================================================================================================================
# INITIAL STATE SPECIFICATION
# ======================================================================================================================================================================================================

DIS_tully1D_1=1; AIS_tully1D_1=1
DIS_tully1D_2=1; AIS_tully1D_2=1
DIS_tully1D_3=0; AIS_tully1D_3=1

DIS_doubleState1D_1=1; AIS_doubleState1D_1=1
DIS_doubleState1D_2=1; AIS_doubleState1D_2=1

DIS_tripleState1D_1=2; AIS_tripleState1D_1=2
DIS_tripleState1D_2=2; AIS_tripleState1D_2=2
DIS_tripleState1D_3=1; AIS_tripleState1D_3=2

# ======================================================================================================================================================================================================
# PROVIDED ARGUMENTS PARSING
# ======================================================================================================================================================================================================

while getopts "c:p:s:e:t:r:l:m:h" ARG; do
    case $ARG in
        c) CORES=$OPTARG;;
        s) PSTEP=$OPTARG;;
        t) TRAJS=$OPTARG;;
        m) MODELS=($OPTARG);;
        h) echo "USAGE: $0 [-c CORES] [-s PSTEP] [-t TRAJS] [-m MODELS] [-h]"; exit 0;;
        *) echo "INVALID OPTION: $ARG"; exit 1;;
    esac
done

echo "CORES: $CORES, PSTEP: $PSTEP, TRAJS: $TRAJS, MODELS: ${MODELS[@]}"

# ======================================================================================================================================================================================================
# EXECUTABLE PATHS
# ======================================================================================================================================================================================================

ACORN=./../../bin/acorn; PLOT_1D=./../../script/plot/plot-1d.py

# ======================================================================================================================================================================================================
# START TEMPLATES
# ======================================================================================================================================================================================================

read -r -d '' TEMPLATE_POTENTIAL_DIABATIC <<- EOM
{
    "adiabatic" : false,
    "limits" : [-16, 16],
    "points" : 512
}
EOM

read -r -d '' TEMPLATE_POTENTIAL_ADIABATIC <<- EOM
{
    "adiabatic" : true,
    "limits" : [-16, 16],
    "points" : 512
}
EOM

read -r -d '' TEMPLATE_EXACT_DYN <<- EOM
{
    "adiabatic" : true,
    "imaginary" : false,
    "iterations" : 500,
    "time_step" : 10,
    "grid" : {
        "limits" : [-48, 144],
        "points" : 4096
    },
    "initial_conditions" : {
        "mass" : 2000,
        "position" : [-15]
    },
    "log_intervals" : {
        "iteration" : 50
    },
    "write" : {}
}
EOM

read -r -d '' TEMPLATE_SH_DYN <<- EOM
{
    "adiabatic" : true,
    "derivative_step" : 0.001,
    "iterations" : 5000,
    "seed" : 1,
    "time_step" : 1,
    "trajectories" : $TRAJS,
    "initial_conditions" : {
        "mass" : 2000,
        "momentum_std" : [1],
        "position_mean" : [-15],
        "position_std" : [0.5]
    },
    "log_intervals" : {
        "iteration" : 500,
        "trajectory" : 100
    },
    "write" : {}
}
EOM

# ======================================================================================================================================================================================================
# END TEMPLATES
# ======================================================================================================================================================================================================

# remove all the previous files
[[ $CLEAN -eq 1 ]] && rm -f *.json *.mat *.png

# generate momenta
(( $(echo "$PSTEP <= 0" | bc -l) )) && MOMENTA=($PSTART) || MOMENTA=($(seq $PSTART $PSTEP $PEND))

# loop over the models
for MODEL in ${MODELS[@]}; do

    # create the json files for potentials
    echo "{}" > potd_${MODEL}.json && echo "{}" > pota_${MODEL}.json

    # fill the json files for potentials
    jq --arg options "${TEMPLATE_POTENTIAL_DIABATIC}"  '. |= . + {"model_potential" : ($options | fromjson)}' "potd_${MODEL}.json" > temp.json && mv temp.json "potd_${MODEL}.json"
    jq --arg options "${TEMPLATE_POTENTIAL_ADIABATIC}" '. |= . + {"model_potential" : ($options | fromjson)}' "pota_${MODEL}.json" > temp.json && mv temp.json "pota_${MODEL}.json"

    # fill the potential in the json files for potentials
    jq --arg potential "${MODEL}" '.model_potential |= . + {"potential" : $potential}' "potd_${MODEL}.json" > temp.json && mv temp.json "potd_${MODEL}.json"
    jq --arg potential "${MODEL}" '.model_potential |= . + {"potential" : $potential}' "pota_${MODEL}.json" > temp.json && mv temp.json "pota_${MODEL}.json"

    # fill the output names in the json files for potentials
    jq --arg potential  "POTENTIAL_DIABATIC_${MODEL}.mat" '.model_potential |= . + {"output" : $potential}' "potd_${MODEL}.json" > temp.json && mv temp.json "potd_${MODEL}.json"
    jq --arg potential "POTENTIAL_ADIABATIC_${MODEL}.mat" '.model_potential |= . + {"output" : $potential}' "pota_${MODEL}.json" > temp.json && mv temp.json "pota_${MODEL}.json"

    # loop over the momentum
    for MOMENTUM in ${MOMENTA[@]}; do

        # create the json files
        echo "{}" > exact_${MODEL}_P=${MOMENTUM}.json && echo "{}" > fssh_${MODEL}_P=${MOMENTUM}.json && echo "{}" > kfssh_${MODEL}_P=${MOMENTUM}.json && echo "{}" > lzsh_${MODEL}_P=${MOMENTUM}.json

        # extract the initial states
        DIS=$(eval echo \$"DIS_${MODEL}"); AIS=$(eval echo \$"AIS_${MODEL}")

        # fill the json files with the templated dynamics
        jq --arg dynamics "${TEMPLATE_EXACT_DYN}" '. |= . + {"quantum_dynamics"   : ($dynamics | fromjson)}' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg dynamics "${TEMPLATE_SH_DYN}"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg dynamics "${TEMPLATE_SH_DYN}"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg dynamics "${TEMPLATE_SH_DYN}"    '. |= . + {"classical_dynamics" : ($dynamics | fromjson)}'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the potential
        jq --arg potential "${MODEL}" '.quantum_dynamics   |= . + {"potential" : $potential}' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg potential "${MODEL}" '.classical_dynamics |= . + {"potential" : $potential}'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg potential "${MODEL}" '.classical_dynamics |= . + {"potential" : $potential}' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg potential "${MODEL}" '.classical_dynamics |= . + {"potential" : $potential}'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the initial state
        jq --arg state "${DIS}" '.quantum_dynamics  .initial_conditions.state |= ($state | tonumber)' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg state "${AIS}" '.classical_dynamics.initial_conditions.state |= ($state | tonumber)'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg state "${AIS}" '.classical_dynamics.initial_conditions.state |= ($state | tonumber)' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg state "${AIS}" '.classical_dynamics.initial_conditions.state |= ($state | tonumber)'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the momentum
        jq --arg momentum "${MOMENTUM}" '.quantum_dynamics  .initial_conditions.momentum      |= [($momentum | tonumber)]' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg momentum "${MOMENTUM}" '.classical_dynamics.initial_conditions.momentum_mean |= [($momentum | tonumber)]'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg momentum "${MOMENTUM}" '.classical_dynamics.initial_conditions.momentum_mean |= [($momentum | tonumber)]' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg momentum "${MOMENTUM}" '.classical_dynamics.initial_conditions.momentum_mean |= [($momentum | tonumber)]'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json with classical dynamics type
        jq '.classical_dynamics.type |=  "fssh"'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq '.classical_dynamics.type |= "kfssh"' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq '.classical_dynamics.type |=  "lzsh"'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the population output
        jq --arg path "POPULATION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      '.quantum_dynamics  .write.population      |= $path' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  '.classical_dynamics.write.population_mean |= $path'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" '.classical_dynamics.write.population_mean |= $path' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POPULATION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  '.classical_dynamics.write.population_mean |= $path'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the kinetic energy output
        jq --arg path "KINETIC_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat"      '.quantum_dynamics  .write.kinetic_energy      |= $path' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  '.classical_dynamics.write.kinetic_energy_mean |= $path'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" '.classical_dynamics.write.kinetic_energy_mean |= $path' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "KINETIC_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  '.classical_dynamics.write.kinetic_energy_mean |= $path'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the potential energy output
        jq --arg path "POTENTIAL_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat"      '.quantum_dynamics  .write.potential_energy      |= $path' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  '.classical_dynamics.write.potential_energy_mean |= $path'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" '.classical_dynamics.write.potential_energy_mean |= $path' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POTENTIAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  '.classical_dynamics.write.potential_energy_mean |= $path'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the total energy output
        jq --arg path "TOTAL_ENERGY_${MODEL}_P=${MOMENTUM}_EXACT.mat"      '.quantum_dynamics  .write.total_energy      |= $path' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  '.classical_dynamics.write.total_energy_mean |= $path'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" '.classical_dynamics.write.total_energy_mean |= $path' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "TOTAL_ENERGY_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  '.classical_dynamics.write.total_energy_mean |= $path'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the position output
        jq --arg path "POSITION_${MODEL}_P=${MOMENTUM}_EXACT.mat"      '.quantum_dynamics  .write.position      |= $path' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  '.classical_dynamics.write.position_mean |= $path'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" '.classical_dynamics.write.position_mean |= $path' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "POSITION_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  '.classical_dynamics.write.position_mean |= $path'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the momentum output
        jq --arg path "MOMENTUM_${MODEL}_P=${MOMENTUM}_EXACT.mat"      '.quantum_dynamics  .write.momentum      |= $path' "exact_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "exact_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  '.classical_dynamics.write.momentum_mean |= $path'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" '.classical_dynamics.write.momentum_mean |= $path' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "MOMENTUM_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  '.classical_dynamics.write.momentum_mean |= $path'  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"

        # fill the json files with the FSSH coefficient output
        jq --arg path "FSSH_COEFFICIENT_MEAN_${MODEL}_P=${MOMENTUM}_FSSH.mat"  '.classical_dynamics.write.fssh_coefficient_mean |= $path'  "fssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "fssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "FSSH_COEFFICIENT_MEAN_${MODEL}_P=${MOMENTUM}_KFSSH.mat" '.classical_dynamics.write.fssh_coefficient_mean |= $path' "kfssh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json "kfssh_${MODEL}_P=${MOMENTUM}.json"
        jq --arg path "FSSH_COEFFICIENT_MEAN_${MODEL}_P=${MOMENTUM}_LZSH.mat"  '.classical_dynamics.write.fssh_coefficient_mean |= null '  "lzsh_${MODEL}_P=${MOMENTUM}.json" > temp.json && mv temp.json  "lzsh_${MODEL}_P=${MOMENTUM}.json"
    done
done
