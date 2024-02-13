#!/bin/bash

# USAGE: wraporca.sh METHOD SYSTEM BASIS CHARGE MULT
if [ "$#" -ne 5 ]; then
    echo "USAGE: wraporca.sh METHOD SYSTEM BASIS CHARGE MULT"; exit 1
fi

# create a folder
mkdir -p .orca

# specify the input file
cat << EOT > .orca/orca.inp
! ${1^^} ${3^^} HCORE NOFROZENCORE

*xyzfile $4 $5 ../$2
EOT

# run the calculation
cd .orca && orca orca.inp; cd .. && rm -rf .orca
