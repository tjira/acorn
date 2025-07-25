#!/bin/bash

usage() {
  cat <<EOF
Usage: $(basename $0) [options]

Options:
  -a <active> Number of active orbitals. (default: ${NACTIVE})
  -b <basis>  Basis set used to specify the atomic orbitals. (default: ${BASIS})
  -c <closed> Number of closed orbitals. (default: ${NCLOSED})
  -g <charge> Charge of the system. (default: ${CHARGE})
  -k <nstate> Number of states to compute. (default: ${NSTATE})
  -m <method> Method to perform. (default: ${METHOD})
  -p <nspin>  The number associated with the spin states: 0 for singlet, 1 for doublet, 2 for triplet, etc. (default: ${NSPIN})
  -r <retain> Retain the .bagel directory and create .all files. (default: ${RETAIN})
  -s <system> System file. (default: ${METHOD})
  -h          Display this help message and exit.
EOF

}; SYSTEM="molecule.xyz"; BASIS="svp"; CHARGE=0; NSPIN=0; METHOD="casscf"; NSTATE=2; NACTIVE=2; NCLOSED=4; RETAIN=0

while getopts "a:b:c:g:k:m:p:rs:h" OPT; do case "$OPT" in
  a ) NACTIVE="$OPTARG" ;;
  b ) BASIS="$OPTARG" ;;
  c ) NCLOSED="$OPTARG" ;;
  g ) CHARGE="$OPTARG" ;;
  k ) NSTATE="$OPTARG" ;;
  m ) METHOD="$OPTARG" ;;
  p ) NSPIN="$OPTARG" ;;
  r ) RETAIN=1 ;;
  s ) SYSTEM="$OPTARG" ;;
  h ) usage && exit 0 ;;
  \? ) usage && exit 1 ;;
esac done

mkdir -p .bagel; METHOD=${METHOD^^}; NATOM=$(awk 'END {print NR - 2}' $SYSTEM)

cat << EOT > .bagel/bagel.json
{ "bagel" : [

{
  "title" : "molecule",
  "basis" : "$BASIS",
  "df_basis" : "$BASIS-jkfit",
  "angstrom" : true,
  "geometry" : [
EOT

awk 'NR == 1 {A = $1} NR > 2 {printf("{\"atom\" : \"%s\", \"xyz\" : [%20.14f, %20.14f, %20.14f]}%s\n", $1, $2, $3, $4, NR == A + 2 ? "" : ",")}' $SYSTEM >> .bagel/bagel.json

cat << EOT >> .bagel/bagel.json
  ]
},
EOT

[[ -f .bagel/orbitals.archive ]] && cat << EOT >> .bagel/bagel.json
{
  "title" : "load_ref",
  "file" : "orbitals",
  "continue_geom" : false
},
EOT

cat << EOT >> .bagel/bagel.json
{
  "title" : "forces",
  "grads" : [
EOT

for I in $(seq 0 $((NSTATE - 1))); do
    awk -v I=$I -v NSTATE=$NSTATE 'BEGIN {printf("{\"title\" : \"force\", \"target\" : %d}%s", I, NSTATE > 1 ? "," : "")}' >> .bagel/bagel.json
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
    "nstate" : $NSTATE,
    "maxiter" : 2000
  }]
},
{
  "title" : "save_ref",
  "file" : "orbitals"
}
]}
EOT

jq . .bagel/bagel.json > .bagel/bagel.json.tmp && mv .bagel/bagel.json.tmp .bagel/bagel.json && rm -f ENERGY.mat GRADIENT.mat NACV.mat

cd .bagel && BAGEL bagel.json | tee bagel.out ; cd ..

if grep -q "EXCEPTION RAISED"             .bagel/bagel.out; then exit 1; fi
if grep -q "Max iteration reached in SCF" .bagel/bagel.out; then exit 1; fi

awk -v NATOM=$NATOM -v NSTATE=$NSTATE 'BEGIN {print NSTATE, 1}'                                     >   ENERGY.mat
awk -v NATOM=$NATOM -v NSTATE=$NSTATE 'BEGIN {print 3 * NATOM * NSTATE, 1}'                         > GRADIENT.mat
awk -v NATOM=$NATOM -v NSTATE=$NSTATE 'BEGIN {print 3 * NATOM * (NSTATE * NSTATE - NSTATE) / 2, 1}' >     NACV.mat

awk '{printf("%20.14f\n", $1)}' .bagel/ENERGY.out >> ENERGY.mat

for I in $(seq 0 $((NSTATE - 1))); do
    awk 'NF != 0 && NR > 1 {printf("%20.14f\n%20.14f\n%20.14f\n", $2, $3, $4)}' .bagel/FORCE_$I.out >> GRADIENT.mat
done

for I in $(seq 0 $((NSTATE - 1))); do
    for J in $(seq $((I + 1)) $((NSTATE - 1))); do
        awk 'NF != 0 && NR > 1 {printf("%20.14f\n%20.14f\n%20.14f\n", $2, $3, $4)}' .bagel/NACME_${I}_$J.out >> NACV.mat
    done
done

if [[ $RETAIN -eq 1 ]]; then
    cat $SYSTEM >> $SYSTEM.all && cat .bagel/bagel.out >> bagel.out.all
else
    rm -rf .bagel
fi
