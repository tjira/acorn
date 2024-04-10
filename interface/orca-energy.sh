#!/bin/bash

# USAGE: orca.sh BASIS CHARGE MULT METHOD
if [ "$#" -ne 4 ]; then
    echo "USAGE: orca.sh BASIS CHARGE MULT METHOD"; exit 1
fi

# specify the input file
cat << EOT > orca.inp
! ${4^^} ${1^^}

*xyzfile $2 $3 molecule.xyz
EOT

# run the calculation
orca orca.inp > orca.out

# extract the final energy
grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}' > energy.dat
