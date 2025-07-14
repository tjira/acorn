#!/bin/bash

cat << EOF > input_1D.json
{
    "quantum_dynamics" : {
        "adiabatic" : true,
        "iterations" : 1,
        "mode" : [0, 1],
        "time_step" : 1,
        "grid" : {
            "limits" : [0.001, 32],
            "points" : 2
        },
        "hamiltonian" : {
            "dims" : 1,
            "file" : "U_1D.dat"
        },
        "initial_conditions" : {
            "mass" : 1,
            "momentum" : [0],
            "position" : [15],
            "state" : 0,
            "gamma" : 2,
            "adiabatic" : true
        },
        "log_intervals" : {
            "iteration" : 500
        },
        "write" : {
            "population" : "POPULATION.mat",
            "position" : "POSITION.mat"
        }
    }
}
EOF

MASS=(1713.8065339025395 3209.8064174768374 4523.60441242962); ATOM=("H" "D" "T")

for I in $(seq 0.035 0.001 1); do
    for J in "${!MASS[@]}"; do

        M=${MASS[$J]}; A=${ATOM[$J]}; cp input_1D.json "input_1D_${A}_E=$I.json"

        sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_1D_${A}_E=$I.json"

        sed -i 's/"iterations" : 1/"iterations" : 7000/' "input_1D_${A}_E=$I.json"

        sed -i 's/"time_step" : 1/"time_step" : 1/' "input_1D_${A}_E=$I.json"

        sed -i 's/"points" : 2/"points" : 4096/' "input_1D_${A}_E=$I.json"

        sed -i 's/"momentum" : \[0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"'\]/' "input_1D_${A}_E=$I.json"

        sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"
    done
done

for I in $(seq 1.05 0.05 10); do
    for J in "${!MASS[@]}"; do

        M=${MASS[$J]}; A=${ATOM[$J]}; cp input_1D.json "input_1D_${A}_E=$I.json"

        sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_1D_${A}_E=$I.json"

        sed -i 's/"iterations" : 1/"iterations" : 4500/' "input_1D_${A}_E=$I.json"

        sed -i 's/"time_step" : 1/"time_step" : 1/' "input_1D_${A}_E=$I.json"

        sed -i 's/"points" : 2/"points" : 4096/' "input_1D_${A}_E=$I.json"

        sed -i 's/"momentum" : \[0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"'\]/' "input_1D_${A}_E=$I.json"

        sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"
    done
done

for I in $(seq 11 1 100); do
    for J in "${!MASS[@]}"; do

        M=${MASS[$J]}; A=${ATOM[$J]}; cp input_1D.json "input_1D_${A}_E=$I.json"

        sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_1D_${A}_E=$I.json"

        sed -i 's/"iterations" : 1/"iterations" : 4500/' "input_1D_${A}_E=$I.json"

        sed -i 's/"time_step" : 1/"time_step" : 0.1/' "input_1D_${A}_E=$I.json"

        sed -i 's/"points" : 2/"points" : 16384/' "input_1D_${A}_E=$I.json"

        sed -i 's/"momentum" : \[0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"'\]/' "input_1D_${A}_E=$I.json"

        sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_1D_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"
    done
done

rm input_1D.json
