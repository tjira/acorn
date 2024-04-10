#!/bin/bash

# define the systems, bases, and methods
SYSTEMS=("H2" "HF" "HCl" "water" "ammonia" "methane" "ethylene")
BASES=("sto-3g" "6-31g" "6-311g" "def2-svp" "cc-pvdz")
METHODS=("HF" "MP2")

# print the header
printf "%10s %12s %6s %13s %13s %8s\n" "SYSTEM" "BASE" "METHOD" "ACORN" "ORCA" "DIFFERENCE"

# loop over the systems, bases, and methods
for SYSTEM in ${SYSTEMS[@]}; do
    for BASE in ${BASES[@]}; do
        for METHOD in ${METHODS[@]}; do
            # run the calculations and extract the energy
            ACORN=$(./wrapper/wrapacorn.sh "./example/molecule/$SYSTEM.xyz" "$BASE" 0 1 "$METHOD" | grep "ENERGY" | tail -n 1 | awk '{print $NF}')
            ORCA=$(./wrapper/wraporca.sh "./example/molecule/$SYSTEM.xyz" "$BASE" 0 1 "$METHOD" | grep "ENERGY" | tail -n 1 | awk '{print $NF}')

            # print the results
            printf "%10s %12s %6s %13.8f %13.8f %.4e\n" "$SYSTEM" "$BASE" "$METHOD" "$ACORN" "$ORCA" "$(echo "sqrt(($ACORN - $ORCA)^2)" | bc -l)"
        done
    done
done
