#!/bin/bash

cat << EOF > input_1D.json
{
    "quantum_dynamics" : {
        "adiabatic" : true,
        "iterations" : 4000,
        "mode" : [0, 1],
        "time_step" : 1,
        "grid" : {
            "limits" : [0.001, 20],
            "points" : 2048
        },
        "hamiltonian" : {
            "dims" : 1,
            "file" : "U.dat"
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

for I in 0.11 $(seq 0.2 0.1 10); do
    for J in "${!MASS[@]}"; do

        M=${MASS[$J]}; A=${ATOM[$J]}; cp input_1D.json "input_1D_${A}_E=$I.json"

        sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_1D_${A}_E=$I.json"

        sed -i 's/"momentum" : \[0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"'\]/' "input_1D_${A}_E=$I.json"

        sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"
    done
done

rm input_1D.json
