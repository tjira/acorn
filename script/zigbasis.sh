#!/bin/bash

BASES=(
    "sto-3g"
    "6-31g"
    "cc-pvdz"
)

for BASIS in ${BASES[@]}; do
    bse get-basis ${BASIS} json --unc-spdf > src/basis/${BASIS}.json
done
