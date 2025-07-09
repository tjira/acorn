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

for I in 0.035 $(seq 0.04 0.035 0.11); do
    for J in "${!MASS[@]}"; do

        M=${MASS[$J]}; A=${ATOM[$J]}; cp input_1D.json "input_1D_${A}_E=$I.json"

        sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_1D_${A}_E=$I.json"

        sed -i 's/"iterations" : 1/"iterations" : 7000/' "input_1D_${A}_E=$I.json"

        sed -i 's/"time_step" : 1/"time_step" : 1/' "input_1D_${A}_E=$I.json"

        sed -i 's/"points" : 2/"points" : 4096/' "input_1D_${A}_E=$I.json"

        sed -i 's/"momentum" : \[0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"'\]/' "input_1D_${A}_E=$I.json"

        sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"
    done
done

for I in 0.2 $(seq 0.5 1.5 11); do
    for J in "${!MASS[@]}"; do

        M=${MASS[$J]}; A=${ATOM[$J]}; cp input_1D.json "input_1D_${A}_E=$I.json"

        sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_1D_${A}_E=$I.json"

        sed -i 's/"iterations" : 1/"iterations" : 4000/' "input_1D_${A}_E=$I.json"

        sed -i 's/"time_step" : 1/"time_step" : 1/' "input_1D_${A}_E=$I.json"

        sed -i 's/"points" : 2/"points" : 4096/' "input_1D_${A}_E=$I.json"

        sed -i 's/"momentum" : \[0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"'\]/' "input_1D_${A}_E=$I.json"

        sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"
    done
done

for I in 15 $(seq 20 10 100); do
    for J in "${!MASS[@]}"; do

        M=${MASS[$J]}; A=${ATOM[$J]}; cp input_1D.json "input_1D_${A}_E=$I.json"

        sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_1D_${A}_E=$I.json"

        sed -i 's/"iterations" : 1/"iterations" : 4000/' "input_1D_${A}_E=$I.json"

        sed -i 's/"time_step" : 1/"time_step" : 0.1/' "input_1D_${A}_E=$I.json"

        sed -i 's/"points" : 2/"points" : 16384/' "input_1D_${A}_E=$I.json"

        sed -i 's/"momentum" : \[0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"'\]/' "input_1D_${A}_E=$I.json"

        sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"

        sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_E='"$I"'.mat"/' "input_1D_${A}_E=$I.json"
    done
done

rm input_1D.json
