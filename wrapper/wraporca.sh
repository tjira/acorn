#!/bin/bash

# USAGE: wraporca.sh SYSTEM BASIS CHARGE MULT METHOD
if [ "$#" -ne 5 ]; then
    echo "USAGE: wraporca.sh SYSTEM BASIS CHARGE MULT METHOD"; exit 1
fi

# create a folder
mkdir -p .orca

# specify the input file
cat << EOT > .orca/orca.inp
! ${5^^} ${2^^} HCORE KDIIS NOFROZENCORE

*xyzfile $3 $4 ../$1
EOT

# run the calculation
cd .orca && orca orca.inp; cd .. && rm -rf .orca
