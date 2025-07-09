#!/bin/bash

export OMP_NUM_THREADS=1; ./input.sh && PATH="$PWD:$PATH" acorn *.json && ./plot.sh
