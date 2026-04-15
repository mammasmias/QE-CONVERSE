#!/bin/bash

# clean directory
rm -rf *.out scratch

$PW -in pw_scf.in > scf.out 2>&1

mpirun --oversubscribe -np 4 $QECONVERSE -nk 2 < gtensor_1.in > gtensor_1.out
mpirun --oversubscribe -np 4 $QECONVERSE -nk 2 < gtensor_2.in > gtensor_2.out

python3 ../check_gtensor.py --files gtensor_1.out gtensor_2.out --refdir reference
