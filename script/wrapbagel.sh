#!/bin/bash

usage() {
  cat <<EOF
Usage: $(basename $0) [options]

Options:
  -s <system> System file. (default: ${METHOD})
  -b <basis>  Basis set used to specify the atomic orbitals. (default: ${BASIS})
  -c <charge> Charge of the system. (default: ${CHARGE})
  -p <nspin>  The number associated with the spin states: 0 for singlet, 1 for doublet, 2 for triplet, etc. (default: ${NSPIN})
  -m <method> Method to perform. (default: ${METHOD})
  -h          Display this help message and exit.
EOF

}; SYSTEM="molecule.xyz"; BASIS="svp"; CHARGE=0; NSPIN=0; METHOD="sa2-casscf(2,4)"

while getopts "b:c:p:m:s:h" OPT; do case "$OPT" in
  b ) BASIS="$OPTARG" ;;
  c ) CHARGE="$OPTARG" ;;
  m ) METHOD="$OPTARG" ;;
  p ) NSPIN="$OPTARG" ;;
  s ) SYSTEM="$OPTARG" ;;
  h ) usage && exit 0 ;;
  \? ) usage && exit 1 ;;
esac done

mkdir -p .bagel; METHOD=${METHOD^^}; NATOM=$(awk 'END {print NR - 2}' $SYSTEM); NSTATE=$(sed -En 's/^SA([0-9]+)-CASSCF.*$/\1/p' <<< $METHOD)

NACTIVE=$(sed -En 's/^SA[0-9]+-CASSCF\(([0-9]+),.*$/\1/p' <<< $METHOD); NCLOSED=$(sed -En 's/^SA[0-9]+-CASSCF\([0-9]+,([0-9]+)\).*$/\1/p' <<< $METHOD)

cat << EOT > .bagel/bagel.json
{ "bagel" : [

{
  "title" : "molecule",
  "basis" : "$BASIS",
  "df_basis" : "$BASIS-jkfit",
  "angstrom" : true,
  "geometry" : [
EOT

awk 'NR == 1 {A = $1} NR > 2 {printf("{\"atom\" : \"%s\", \"xyz\" : [%20.14f, %20.14f, %20.14f]}%s\n", $1, $2, $3, $4, NR == A + 2 ? "" : ",")}' molecule.xyz >> .bagel/bagel.json

cat << EOT >> .bagel/bagel.json
  ]
},
EOT

cat << EOT >> .bagel/bagel.json
{
  "title" : "forces",
  "grads" : [
EOT

for I in $(seq 0 $((NSTATE - 1))); do
    awk -v I=$I 'BEGIN {printf("{\"title\" : \"force\", \"target\" : %d},", I)}' >> .bagel/bagel.json
done

for I in $(seq 0 $((NSTATE - 1))); do
    for J in $(seq $((I + 1)) $((NSTATE - 1))); do
        awk -v I=$I -v J=$J -v NSTATE=$NSTATE 'BEGIN {printf("{\"title\" : \"nacme\", \"target\" : %d, \"target2\" : %d}%s", I, J, I == NSTATE - 2 && J == NSTATE - 1 ? "" : ",")}' >> .bagel/bagel.json
    done
done

cat << EOT >> .bagel/bagel.json
  ],
  "export" : true,
  "method" : [{
    "title" : "casscf",
    "charge" : $CHARGE,
    "nspin" : $NSPIN,
    "nact" : $NACTIVE,
    "nclosed" : $NCLOSED,
    "nstate" : $NSTATE
  }]
}
]}
EOT

jq . .bagel/bagel.json > .bagel/bagel.json.tmp && mv .bagel/bagel.json.tmp .bagel/bagel.json

cd .bagel && BAGEL bagel.json ; cd ..

awk -v NATOM=$NATOM -v NSTATE=$NSTATE 'BEGIN {print NSTATE + 3 * NATOM * NSTATE + 3 * NATOM * (NSTATE * NSTATE - NSTATE) / 2}' > ENGRADNACV.mat

awk '{printf("%20.14f\n", $1)}' .bagel/ENERGY.out >> ENGRADNACV.mat

for I in $(seq 0 $((NSTATE - 1))); do
    awk 'NF != 0 && NR > 1 {printf("%20.14f\n%20.14f\n%20.14f\n", $2, $3, $4)}' .bagel/FORCE_$I.out >> ENGRADNACV.mat
done

for I in $(seq 0 $((NSTATE - 1))); do
    for J in $(seq $((I + 1)) $((NSTATE - 1))); do
        awk 'NF != 0 && NR > 1 {printf("%20.14f\n%20.14f\n%20.14f\n", $2, $3, $4)}' .bagel/NACME_${I}_$J.out >> ENGRADNACV.mat
    done
done

rm -rf .bagel
