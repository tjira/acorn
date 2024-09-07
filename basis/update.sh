#!/bin/bash

# define the bases to download
BASES=(
    "sto-2g"
    "sto-3g"
    "sto-4g"
    "sto-5g"
    "sto-6g"
    "6-31g"
    "6-31g*"
    "6-31g**"
    "6-31+g"
    "6-31+g*"
    "6-31+g**"
    "6-31++g"
    "6-31++g*"
    "6-31++g**"
    "6-311g"
    "6-311g*"
    "6-311g**"
    "6-311+g"
    "6-311+g*"
    "6-311+g**"
    "6-311++g"
    "6-311++g*"
    "6-311++g**"
    "def2-svp"
    "def2-tzvp"
    "def2-qzvp"
    "def2-svpd"
    "def2-tzvpd"
    "def2-qzvpd"
    "def2-tzvpp"
    "def2-qzvpp"
    "def2-tzvppd"
    "def2-qzvppd"
    "cc-pvdz"
    "cc-pvtz"
    "cc-pvqz"
    "aug-cc-pvdz"
    "aug-cc-pvtz"
    "aug-cc-pvqz"
)

# download the basis sets
for BASE in "${BASES[@]}"; do
    bse -o "$(echo "${BASE}" | sed "s/*/s/g ; s/+/p/g ; s/aug-/augmented-/g").g94" get-basis "${BASE}" gaussian94 --unc-spdf
done
