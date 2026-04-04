#!/usr/bin/env bash

# clean directory
rm -rf *.out scratch

mpirun -np 8 ~/src/q-e/bin/pw.x < pw.scf.in > scf.out

mpirun -np 8 ~/src/QE-CONVERSE/bin/qe-converse.x < li1x.in > li1x.out
mpirun -np 8 ~/src/QE-CONVERSE/bin/qe-converse.x < li1y.in > li1y.out
mpirun -np 8 ~/src/QE-CONVERSE/bin/qe-converse.x < li1z.in > li1z.out
