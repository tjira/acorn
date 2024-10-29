#!/bin/bash

# check the number of arguments
[ $# -ne 1 ] && echo "USAGE: $0 MODE" && exit 1

# get the latest tag
TAG=$(git describe --tags --abbrev=0)

# download the libraries
rm -rf external && mkdir -p external && wget -O "external/libs-${1,,}.tar.gz" "https://github.com/tjira/acorn/releases/download/$TAG/acorn-libs_linux-${1,,}_x86-64.tar.gz"

# extract the libraries
cd external && tar -xzvf "libs-${1,,}.tar.gz" && cd ..
