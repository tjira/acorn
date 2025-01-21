#!/bin/bash

BIN=../../zig-out/x86_64-linux/example; PLOT=../../python/plot.py; export PLOT;

zig build --release=fast && rm -rf result *.mat *.png && mkdir -p result && for FILE in *.zig; do
    [[ "$#" -eq 1 ]] && [[ $1 != $FILE ]] && continue; $BIN/${FILE%.*} && [[ -f ./analyze/${FILE%.*}.sh ]] && ./analyze/${FILE%.*}.sh && mkdir result/${FILE%.*} &&  mv *.mat *.png result/${FILE%.*};
done
