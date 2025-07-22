#!/bin/bash

usage() {
  cat <<EOF
Usage: $(basename $0) [options]

Options:
  -s <system>       System file. (default: ${SYSTEM})
  -b <basis>        Basis set used to specify the atomic orbitals. (default: ${BASIS})
  -c <charge>       Charge of the system. (default: ${CHARGE})
  -p <multiplicity> Spin multiplicity of the system. (default: ${MULTIPLICITY})
  -m <method>       Method to perform. (default: ${METHOD})
  -h                Display this help message and exit.
EOF

}; SYSTEM="molecule.xyz"; BASIS="sto-3g"; CHARGE=0; MULTIPLICITY=1; METHOD="hf"

while getopts "b:c:p:m:s:h" OPT; do case "$OPT" in
  b ) BASIS="$OPTARG" ;;
  c ) CHARGE="$OPTARG" ;;
  m ) METHOD="$OPTARG" ;;
  p ) MULTIPLICITY="$OPTARG" ;;
  s ) SYSTEM="$OPTARG" ;;
  h ) usage && exit 0 ;;
  \? ) usage && exit 1 ;;
esac done

mkdir -p .orca; METHOD=${METHOD^^}; NATOM=$(awk 'END {print NR - 2}' $SYSTEM)

if [[ $METHOD == "FCI" ]]; then
    METHOD=""; OPTIONS="%CASSCF DOFCI TRUE; END"
fi

cat << EOT > .orca/orca.inp
! $METHOD ${BASIS^^} HCORE KDIIS NOFROZENCORE TIGHTSCF ENGRAD LARGEPRINT

$OPTIONS

*xyzfile $CHARGE $MULTIPLICITY ../$SYSTEM
EOT

sed -i 's/   /  /g ; s/  / /g' .orca/orca.inp

cd .orca && orca orca.inp ; cd ..

if [[ -f .orca/orca.engrad ]]; then
    awk -v NATOM=$NATOM 'BEGIN {print 3*NATOM+1} NR==8 || (NR>11 && NR<=3*NATOM+11) {printf "%20.14f\n", $1}' .orca/orca.engrad > ENGRAD.mat
fi

rm -rf .orca
