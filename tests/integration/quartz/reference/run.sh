#!/bin/bash

# clean directory
rm -rf *.out scratch

mpirun -np 4 pw.x -in pw_scf.in > scf.out 2>&1

mpirun -np 4 ../../../../src/qe-converse.x < Si1x.in > Si1x.out
mpirun -np 4 ../../../../src/qe-converse.x < O4z.in > O4z.out

