#!/bin/bash

# create folders
mkdir -p external && mkdir -p external/include

# download argparse
wget -q -O external/include/argparse.hpp https://raw.githubusercontent.com/p-ranav/argparse/master/include/argparse/argparse.hpp

# download json
wget -q -O external/include/json.hpp https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp

# download exprtk
wget -q -O external/include/exprtk.hpp https://raw.githubusercontent.com/ArashPartow/exprtk/master/exprtk.hpp
