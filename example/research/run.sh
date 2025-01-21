#!/bin/bash

BIN=../../zig-out/x86_64-linux/example; PLOT=../../python/plot.py; export PLOT;

rm -rf result *.mat *.png && mkdir -p result

for FILE in *.zig; do
    $BIN/${FILE%.*} && [[ -f ./analyze/${FILE%.*}.sh ]] && ./analyze/${FILE%.*}.sh && mkdir result/${FILE%.*} &&  mv *.mat *.png result/${FILE%.*};
done
