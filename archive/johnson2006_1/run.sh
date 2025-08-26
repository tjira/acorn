#!/bin/bash

export OMP_NUM_THREADS=1; PATH="$PWD:$PATH" acorn *.json
