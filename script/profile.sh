#!/bin/bash

zig build

valgrind --callgrind-out-file=callgrind.out --tool=callgrind acorn $1

gprof2dot --format=callgrind --output=profile.dot callgrind.out

dot -T pdf profile.dot -o profile.pdf

rm -f callgrind.out profile.dot
