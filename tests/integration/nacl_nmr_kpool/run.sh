#!/usr/bin/env bash

# clean directory
rm -rf *.out scratch

mpirun --oversubscribe -np 4 $PW -in pw.scf.in > scf.out 2>&1

mpirun --oversubscribe -np 8 $QECONVERSE -nk 4 < Na1y.in > Na1y.out
mpirun --oversubscribe -np 8 $QECONVERSE -nk 4 < Cl2z.in > Cl2z.out
wait

python3 ../check_nmr.py \
    --files Na1y.out Cl2z.out \
    --refdir reference
