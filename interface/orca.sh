#!/bin/bash

# create the input file
cat << EOT > orca.inp
! ${4^^} ENGRAD ${3^^}

*xyzfile $1 $2 molecule.xyz
EOT

# numerical gradients for some methods
if [[ ${4^^} == "CISD" ]] || [[ ${4^^} == "CISD(T)" ]] || [[ ${4^^} == "CCSD" ]] || [[ ${4^^} == "CCSD(T)" ]]; then
    sed -i 's/ENGRAD/ENGRAD NUMGRAD/' orca.inp
fi

# run the calculation
orca orca.inp > orca.out

# extract the final energy
grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}' > energy.dat

# extract the gradient
sed -n "12,$(echo "11 + $(head -n 1 molecule.xyz) * 3" | bc -l)p" orca.engrad > gradient.dat
