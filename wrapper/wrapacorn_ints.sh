#!/bin/bash

# USAGE: wrapacorn.sh SYSTEM BASIS
if [ "$#" -ne 2 ]; then
    echo "USAGE: wrapacorn.sh SYSTEM BASIS"; exit 1
fi

# create the input file
jq -n --arg file "$1" --arg basis "$2" '{molecule: $ARGS.named, integral: {export: {overlap: true, kinetic: true, nuclear: true, coulomb: true}}}' > input.json

# run the calculation and clean up
acorn input.json; rm input.json
