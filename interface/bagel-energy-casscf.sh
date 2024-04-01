#!/bin/bash

# USAGE: bagel.sh CHARGE MULT BASIS
if [ "$#" -ne 5 ]; then
    echo "USAGE: bagel.sh BASIS CHARGE MULT NSTATE STATE"; exit 1
fi

# specify the molecule in the input file
cat << EOT > bagel.json
{ "bagel" : [

{
    "title" : "molecule",
    "basis" : "$1",
    "df_basis" : "$1-jkfit",
    "angstrom" : true,
    "geometry" : [
EOT

# add the coordinates
while read line; do
    echo "        { \"atom\" : \"$(echo $line | awk '{print $1}')\", \"xyz\" : [$(echo $line | awk '{printf "%20.14f, %20.14f, %20.14f", $2, $3, $4}')] }," >> bagel.json
done <<< $(tail -n +3 molecule.xyz)

# remove the last comma
truncate -s -2 bagel.json && echo "" >> bagel.json

# specify the method
cat << EOT >> bagel.json
    ]
},

{
  "title" : "forces",
  "grads" : [
    { "title" : "force", "target" : $5 }
  ],
  "export" : true,
  "method" : [ {
    "title" : "casscf",
    "nstate" : $4,
    "nclosed" : 4,
    "nact" : 2
  } ]
}

]}
EOT

# run the calculation
BAGEL bagel.json > bagel.out

# extract the final energy
mv ENERGY.out energy.dat

# extract the gradient
tail -n +2 "FORCE_$5.out" | awk 'NF {print $2, $3, $4}' | sed 's/ /\n/g' > gradient.dat
