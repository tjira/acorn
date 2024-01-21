#!/bin/bash

# create the input file
cat << EOT > orca.inp
! $4 ENGRAD $3

*xyzfile $1 $2 molecule.xyz
EOT

# run the calculation
orca orca.inp > orca.out

# extract the final energy
grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}' > energy.dat

# extract the gradient
sed -n "12,$(echo "11 + $(head -n 1 molecule.xyz) * 3" | bc -l)p" orca.engrad > gradient.dat
