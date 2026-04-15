#!/bin/bash

# clean directory
rm -rf *.out scratch

pw.x -in pw_scf.in > scratch/scf.out 2>&1

../../../../src/qe-converse.x < gtensor_1.in > gtensor_1.out
../../../../src/qe-converse.x < gtensor_2.in > gtensor_2.out

