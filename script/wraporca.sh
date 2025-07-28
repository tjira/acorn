#!/bin/bash

usage() {
  cat <<EOF
Usage: $(basename $0) [options]

Options:
  -b <basis>        Basis set used to specify the atomic orbitals. (default: ${BASIS})
  -c <charge>       Charge of the system. (default: ${CHARGE})
  -k <nstate>       Number of states to compute. (default: ${NSTATE})
  -m <method>       Method to perform. (default: ${METHOD})
  -p <multiplicity> Spin multiplicity of the system. (default: ${MULTIPLICITY})
  -r <retain>       Retain the .orca directory and create .all files. (default: ${RETAIN})
  -s <system>       System file. (default: ${SYSTEM})
  -h                Display this help message and exit.
EOF

}; SYSTEM="molecule.xyz"; BASIS="sto-3g"; CHARGE=0; MULTIPLICITY=1; METHOD="hf"; NSTATE=1; RETAIN=0

while getopts "b:c:k:m:p:rs:h" OPT; do case "$OPT" in
  b ) BASIS="$OPTARG" ;;
  c ) CHARGE="$OPTARG" ;;
  k ) NSTATE="$OPTARG" ;;
  m ) METHOD="$OPTARG" ;;
  p ) MULTIPLICITY="$OPTARG" ;;
  r ) RETAIN=1 ;;
  s ) SYSTEM="$OPTARG" ;;
  h ) usage && exit 0 ;;
  \? ) usage && exit 1 ;;
esac done

mkdir -p .orca; METHOD=${METHOD^^}; NATOM=$(awk 'END {print NR - 2}' $SYSTEM); GUESS=$([ -f .orca/orca.gbw ] && echo "MOREAD" || echo "HCORE")

if [[ $METHOD == "FCI" ]]; then
    METHOD=""; OPTIONS="%CASSCF DOFCI TRUE; END"
fi

if [[ -f .orca/orca.gbw ]]; then
    mv .orca/orca.gbw .orca/guess.gbw
fi

cat << EOT > .orca/orca.inp
! $METHOD ${BASIS^^} $GUESS KDIIS NOFROZENCORE TIGHTSCF ENGRAD LARGEPRINT

%moinp "guess.gbw"

*xyzfile $CHARGE $MULTIPLICITY $SYSTEM
EOT

sed -i 's/   /  /g ; s/  / /g' .orca/orca.inp

cd .orca && cp ../$SYSTEM $SYSTEM && orca orca.inp | tee orca.out ; cd .. && rm -f ENERGY.mat GRADIENT.mat NACV.mat

awk -v NATOM=$NATOM 'BEGIN {print 1, 1} NR == 8 {printf "%20.14f\n", $1}' .orca/orca.engrad > ENERGY.mat

awk -v NATOM=$NATOM 'BEGIN {print 3 * NATOM, 1} NR > 11 && NR <= 3 * NATOM + 11 {printf "%20.14f\n", $1}' .orca/orca.engrad > GRADIENT.mat

if [[ $RETAIN -eq 1 ]]; then
    cat $SYSTEM >> $SYSTEM.all && cat .orca/orca.out >> orca.out.all
else
    rm -rf .orca
fi
