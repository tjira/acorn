#!/bin/bash

usage() {
  cat <<EOF
Usage: $(basename $0) [options]

Options:
  -s <system>       System file. (default: ${COUNT})
  -b <basis>        Basis set used to specify the atomic orbitals. (default: ${LOG_INTERVAL})
  -c <charge>       Charge of the system. (default: ${OUTPUT})
  -p <multiplicity> Spin multiplicity of the system. (default: ${START})
  -m <method>       Method to perform. (default: ${START})
  -f                Flag for the hessian and frequency calculation.
  -g                Flag for the gradient calculation.
  -h                Display this help message and exit.
EOF

}; SYSTEM="molecule.xyz"; BASIS="sto-3g"; CHARGE=0; MULTIPLICITY=1; METHOD="hf"; FREQUENCY=0; GRADIENT=0

while getopts "s:b:c:p:m:fgh" OPT; do case "$OPT" in
  s ) SYSTEM="$OPTARG" ;; b ) BASIS="$OPTARG" ;; c ) CHARGE="$OPTARG" ;; p ) MULTIPLICITY="$OPTARG" ;; m ) METHOD="$OPTARG" ;; f ) FREQUENCY=1 ;; g ) GRADIENT=1 ;; h ) usage && exit 0 ;; \? ) usage && exit 1 ;;
esac done

mkdir -p .orca; METHOD=${METHOD^^}; OPTIONS=""

if [[ $METHOD == "FCI" ]]; then
    METHOD=""; OPTIONS="%CASSCF DOFCI TRUE; END"
fi

if [[ $FREQUENCY -eq 1 ]]; then
    FREQUENCY="NUMFREQ"; OPTIONS="%FREQ PROJECTTR FALSE; END"; else FREQUENCY=""
fi

if [[ $GRADIENT -eq 1 ]]; then
    GRADIENT="ENGRAD NUMGRAD"; else GRADIENT=""
fi

cat << EOT > .orca/orca.inp
! $METHOD ${BASIS^^} HCORE KDIIS NOFROZENCORE TIGHTSCF $GRADIENT $FREQUENCY LARGEPRINT

$OPTIONS

*xyzfile $CHARGE $MULTIPLICITY ../$SYSTEM
EOT

sed -i 's/   /  /g ; s/  / /g' .orca/orca.inp

cd .orca && orca orca.inp; cd .. && rm -rf .orca
