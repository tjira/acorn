#!/bin/bash

# USAGE: wraporca.sh SYSTEM BASIS CHARGE MULT METHOD
if [ "$#" -ne 5 ]; then
    echo "USAGE: wraporca.sh SYSTEM BASIS CHARGE MULT METHOD"; exit 1
fi

# create a folder
mkdir -p .orca

# extract the method
METHOD=${5^^}; if [[ ${5^^} == "FCI" ]]; then METHOD=""; fi

# extract the options
OPTIONS=""; if [[ ${5^^} == "FCI" ]]; then OPTIONS="%casscf DoFCI true; end"; fi

# specify the input file
cat << EOT > .orca/orca.inp
! $METHOD ${2^^} HCORE KDIIS NOFROZENCORE

$OPTIONS

*xyzfile $3 $4 ../$1
EOT

# run the calculation
cd .orca && orca orca.inp; cd .. && rm -rf .orca
