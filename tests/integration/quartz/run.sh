#!/bin/bash

# clean directory
rm -rf *.out scratch

mpirun --oversubscribe -np 4 $PW -in pw_scf.in > scf.out 2>&1

mpirun --oversubscribe -np 4 $QECONVERSE < Si1x.in > Si1x.out
mpirun --oversubscribe -np 4 $QECONVERSE < O4z.in > O4z.out

python3 ../check_nmr.py --files Si1x.out O4z.out --refdir reference
