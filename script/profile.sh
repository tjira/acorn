#!/bin/bash

zig build

valgrind --callgrind-out-file=callgrind.out --tool=callgrind acorn $1

gprof2dot --edge-thres=1 --format=callgrind --node-thres=5 --output=profile.dot --root=main.main callgrind.out

dot -T pdf profile.dot -o profile.pdf
