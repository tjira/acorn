#!/bin/bash

BASES=("sto-3g" "6-31g" "6-31g*" "cc-pvdz" "cc-pvtz")

for BASE in "${BASES[@]}"; do
    bse -o "${BASE}.g94" get-basis "${BASE}" gaussian94
done
