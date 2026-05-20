#!/usr/bin/env bash

# clean directory
rm -rf *.out scratch

mpirun --oversubscribe -np 9 pw.x -nk 9 -in pw.scf.in > scf.out 2>&1

mpirun --oversubscribe -np 9 ../../../../src/qe-converse.x -nk 9 < li1x.in > li1x.out

# mpirun --oversubscribe -np 9 ../../../../src/qe-converse.x -nk 9 < li2y.in > li2y.out

mpirun --oversubscribe -np 9 ../../../../src/qe-converse.x -nk 9 < co6z.in > co6z.out


# python3 ../check_nmr.py --files li1x.out li2y.out co6z.out co7x.out --refdir reference
