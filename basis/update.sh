#!/bin/bash

# define the bases to download
BASES=("sto-3g" "6-31g" "6-31g*" "cc-pvdz" "cc-pvtz")

# download the basis sets
for BASE in "${BASES[@]}"; do
    bse -o "${BASE}.g94" get-basis "${BASE}" gaussian94
done
