#!/bin/bash

mpirun -np 4 $PW -in pw_scf.in > scratch/scf.out 2>&1

mpirun -np 4 $QECONVERSE < Si1x.in > Si1x.out
mpirun -np 4 $QECONVERSE < Si1y.in > Si1y.out
mpirun -np 4 $QECONVERSE < Si1z.in > Si1z.out

python3 ../check_nmr.py --files Si1x.out Si1y.out Si1z.out --refdir reference
