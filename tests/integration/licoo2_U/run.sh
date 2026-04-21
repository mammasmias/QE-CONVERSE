#!/usr/bin/env bash

# clean directory
rm -rf *.out scratch

mpirun --oversubscribe -np 9 $PW -nk 9 -in pw.scf.in > scf.out 2>&1

mpirun --oversubscribe -np 9 $QECONVERSE -nk 9 < li1x.in > li1x.out
 
# mpirun --oversubscribe -np 9 $QECONVERSE -nk 9 < li2y.in > li2y.out

# mpirun --oversubscribe -np 9 $QECONVERSE -nk 9 < co6z.in > co6z.out

mpirun --oversubscribe -np 9 $QECONVERSE -nk 9 < co7x.in > co7x.out

python3 ../check_nmr.py --files li1x.out co7x.out --refdir reference
