#!/bin/bash

for F in *.json; do
    echo "export OMP_NUM_THREADS=1; PATH=\"\$PWD:\$PATH\" acorn $F" > "run_$F.sh" && chmod +x "run_$F.sh" && qsub -q kq-64-8 -cwd -V "run_$F.sh"
done
