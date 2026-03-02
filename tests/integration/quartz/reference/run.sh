#!/bin/bash

# clean directory
rm -rf *.out scratch

mpirun -np 4 pw.x -in pw_scf.in > scratch/scf.out 2>&1

mpirun -np 4 ../../../../src/qe-converse.x < Si1x.in > Si1x.out
mpirun -np 4 ../../../../src/qe-converse.x < Si1y.in > Si1y.out
mpirun -np 4 ../../../../src/qe-converse.x < Si1z.in > Si1z.out

