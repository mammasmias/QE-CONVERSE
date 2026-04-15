#!/usr/bin/env bash
# Generate NaCl NMR reference outputs.
# Run this script from the reference/ directory to (re)create reference output files.

# clean directory
rm -rf *.out scratch

mpirun -np 4 pw.x -in pw.scf.in > scf.out 2>&1

# NaCl is a tiny 2-atom system: run qe-converse serial (np=1) so 6 jobs
# can run in parallel without memory-bandwidth contention or MPI overhead.
../../../../bin/qe-converse.x < Na1x.in > Na1x.out &
../../../../bin/qe-converse.x < Na1y.in > Na1y.out &
../../../../bin/qe-converse.x < Na1z.in > Na1z.out &

../../../../bin/qe-converse.x < Cl2x.in > Cl2x.out &
../../../../bin/qe-converse.x < Cl2y.in > Cl2y.out &
../../../../bin/qe-converse.x < Cl2z.in > Cl2z.out &
wait
