#!/bin/bash

zig build -DBENCHMARK -DDEBUG

valgrind --callgrind-out-file=callgrind.out --tool=callgrind $@

gprof2dot --edge-thres=1 --format=callgrind --node-thres=5 --output=profile.dot --root=main callgrind.out

dot -T pdf profile.dot -o profile.pdf
