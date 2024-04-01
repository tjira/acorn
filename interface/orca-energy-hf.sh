#!/bin/bash

# USAGE: orca.sh CHARGE MULT BASIS
if [ "$#" -ne 3 ]; then
    echo "USAGE: orca.sh BASIS CHARGE MULT"; exit 1
fi

# specify the input file
cat << EOT > orca.inp
! HF ${1^^}

*xyzfile $2 $3 molecule.xyz
EOT

# run the calculation
orca orca.inp > orca.out

# extract the final energy
grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}' > energy.dat
