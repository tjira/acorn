#!/bin/bash

# define the bases to download
BASES=("sto-3g" "6-31g" "6-31+g" "6-31g*" "cc-pvdz" "cc-pvtz")

# download the basis sets
for BASE in "${BASES[@]}"; do
    bse -o "$(echo "${BASE}" | sed "s/*/s/g ; s/+/p/g").g94" get-basis "${BASE}" gaussian94 --unc-spdf
done
