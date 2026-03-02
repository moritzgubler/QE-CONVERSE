#!/usr/bin/env bash
# Generate NaCl NMR reference outputs.
# Run this script from the reference/ directory to (re)create reference output files.

# clean directory
rm -rf *.out scratch

mpirun -np 4 pw.x -in pw.scf.in > scf.out 2>&1

mpirun -np 4 ../../../../bin/qe-converse.x < Na1x.in > Na1x.out
mpirun -np 4 ../../../../bin/qe-converse.x < Na1y.in > Na1y.out
mpirun -np 4 ../../../../bin/qe-converse.x < Na1z.in > Na1z.out

mpirun -np 4 ../../../../bin/qe-converse.x < Cl2x.in > Cl2x.out
mpirun -np 4 ../../../../bin/qe-converse.x < Cl2y.in > Cl2y.out
mpirun -np 4 ../../../../bin/qe-converse.x < Cl2z.in > Cl2z.out
