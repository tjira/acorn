#!/bin/bash

cat << EOF > input.py
import numpy as np

U_1D = np.loadtxt("U_1D.dat", skiprows=1)

grid_2D = np.stack([G.ravel() for G in np.meshgrid(*[np.linspace(-U_1D[-1, 0], U_1D[-1, 0], 1024)] * 2, indexing="ij")], axis=-1).reshape(-1, 2)
grid_3D = np.stack([G.ravel() for G in np.meshgrid(*[np.linspace(-U_1D[-1, 0], U_1D[-1, 0],  128)] * 3, indexing="ij")], axis=-1).reshape(-1, 3)

V00_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 1])
V01_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 2])
V10_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 3])
V11_2D = np.interp(np.linalg.norm(grid_2D, axis=1), U_1D[:, 0], U_1D[:, 4])

V00_3D = np.interp(np.linalg.norm(grid_3D, axis=1), U_1D[:, 0], U_1D[:, 1])
V01_3D = np.interp(np.linalg.norm(grid_3D, axis=1), U_1D[:, 0], U_1D[:, 2])
V10_3D = np.interp(np.linalg.norm(grid_3D, axis=1), U_1D[:, 0], U_1D[:, 3])
V11_3D = np.interp(np.linalg.norm(grid_3D, axis=1), U_1D[:, 0], U_1D[:, 4])

np.savetxt("U_2D.mat", np.column_stack((grid_2D, V00_2D, V01_2D, V10_2D, V11_2D)), fmt="%20.14f", header=f"{grid_2D.shape[0]} 6", comments="")
np.savetxt("U_3D.mat", np.column_stack((grid_3D, V00_3D, V01_3D, V10_3D, V11_3D)), fmt="%20.14f", header=f"{grid_3D.shape[0]} 6", comments="")
EOF

python input.py && rm input.py

cat << EOF > input_1D.json
{
    "quantum_dynamics" : {
        "adiabatic" : true,
        "iterations" : 1,
        "mode" : [0, 1],
        "time_step" : 1,
        "grid" : {
            "limits" : [[0.001, 32]],
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

cat << EOF > input_2D.json
{
    "quantum_dynamics" : {
        "adiabatic" : true,
        "iterations" : 1,
        "mode" : [0, 1],
        "time_step" : 5,
        "grid" : {
            "limits" : [[-10, 18], [-10, 10]],
            "points" : 2
        },
        "hamiltonian" : {
            "dims" : 2,
            "file" : "U_2D.mat"
        },
        "initial_conditions" : {
            "mass" : 1,
            "momentum" : [0, 0],
            "position" : [15, 0],
            "state" : 0,
            "gamma" : 2,
            "adiabatic" : true
        },
        "log_intervals" : {
            "iteration" : 100
        },
        "write" : {
            "population" : "POPULATION.mat"
        }
    }
}
EOF

MASS=(1713.8065339025395 3209.8064174768374 4523.60441242962); ATOM=("H" "D" "T")

for I in $(seq 0.035 0.001 0.1); do
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

for I in $(seq 0.11 0.01 1); do
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

for I in $(seq 1.1 0.1 10); do
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

for B in 0; do
    for I in $(seq 0.035 0.005 0.1); do
        for J in "${!MASS[@]}"; do

            M=${MASS[$J]}; A=${ATOM[$J]}; cp input_2D.json "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"iterations" : 1/"iterations" : 1000/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"time_step" : 1/"time_step" : 5/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"points" : 2/"points" : 1024/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"position" : \[15, 0\]/"position" : \[15, '"$B"'\]/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"momentum" : \[0, 0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"', 0\]/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"
        done
    done

    for I in $(seq 0.2 0.1 1); do
        for J in "${!MASS[@]}"; do

            M=${MASS[$J]}; A=${ATOM[$J]}; cp input_2D.json "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"iterations" : 1/"iterations" : 500/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"time_step" : 1/"time_step" : 5/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"points" : 2/"points" : 2048/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"position" : \[15, 0\]/"position" : \[15, '"$B"'\]/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"momentum" : \[0, 0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"', 0\]/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"
        done
    done

    for I in $(seq 2 1 10); do
        for J in "${!MASS[@]}"; do

            M=${MASS[$J]}; A=${ATOM[$J]}; cp input_2D.json "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"iterations" : 1/"iterations" : 500/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"time_step" : 1/"time_step" : 2.5/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"points" : 2/"points" : 4096/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"position" : \[15, 0\]/"position" : \[15, '"$B"'\]/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"momentum" : \[0, 0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"', 0\]/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"
        done
    done

    for I in $(seq 20 10 100); do
        for J in "${!MASS[@]}"; do

            M=${MASS[$J]}; A=${ATOM[$J]}; cp input_2D.json "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"mass" : 1/"mass" : '"$M"'/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"iterations" : 1/"iterations" : 500/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"time_step" : 1/"time_step" : 1/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"points" : 2/"points" : 8192/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"position" : \[15, 0\]/"position" : \[15, '"$B"'\]/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"momentum" : \[0, 0\]/"momentum" : \[-'"$(echo "sqrt(2*$M*$I)" | bc -l)"', 0\]/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"population" : "POPULATION.mat"/"population" : "POPULATION_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"position" : "POSITION.mat"/"position" : "POSITION_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"

            sed -i 's/"density" : "DENSITY.mat"/"density" : "DENSITY_'"$A"'_2D_b='"$B"'_E='"$I"'.mat"/' "input_2D_${A}_b=${B}_E=$I.json"
        done
    done
done

rm input_1D.json input_2D.json
