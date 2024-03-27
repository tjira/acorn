#!/bin/bash

# USAGE: wrapacorn.sh METHOD SYSTEM BASIS CHARGE MULT
if [ "$#" -ne 5 ]; then
    echo "USAGE: wrapacorn.sh SYSTEM BASIS CHARGE MULT METHOD"; exit 1
fi

read -d '' RHF << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/test/test_rhf %s %s)
set_tests_properties(%s PROPERTIES PASS_REGULAR_EXPRESSION "%s")
EOF

read -d '' RMP2 << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/test/test_rmp2 %s %s)
set_tests_properties(%s PROPERTIES PASS_REGULAR_EXPRESSION "%s")
EOF

read -d '' RFCI << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/test/test_rfci %s %s)
set_tests_properties(%s PROPERTIES PASS_REGULAR_EXPRESSION "%s")
EOF

create() {
    # extract the variables from the function arguments
    METHOD="${5^^}"; SYSTEM="$1"; BASIS="${2^^}"; CHARGE="$3"; MULT="$4"

    # create the molecule block
    jq -n --arg file "$SYSTEM" --arg basis "$BASIS" --argjson charge "$CHARGE" --argjson multiplicity "$MULT" '{molecule: $ARGS.named}' > molecule.json

    # create the input file
    if [ "$METHOD" == "RHF" ]; then
        jq '. + {rhf: {}}' molecule.json > input.json
    elif [ "$METHOD" == "UHF" ]; then
        jq '. + {uhf: {}}' molecule.json > input.json
    elif [ "$METHOD" == "RMP2" ]; then
        jq '. + {rhf: {}, rmp: {}}' molecule.json > input.json
    elif [ "$METHOD" == "RFCI" ]; then
        jq '. + {rhf: {}, rci: {}}' molecule.json > input.json
    else
        echo "INVALID METHOD: $METHOD"; exit 1
    fi
}

create "$1" "$2" "$3" "$4" "$5"; acorn input.json; rm molecule.json input.json
