#!/usr/bin/env bash

# clean directory
rm -rf *.out scratch

mpirun -np 4 $PW -in pw.scf.in > scf.out 2>&1

mpirun -np 4 $QECONVERSE < Na1x.in > Na1x.out
mpirun -np 4 $QECONVERSE < Na1y.in > Na1y.out
mpirun -np 4 $QECONVERSE < Na1z.in > Na1z.out

mpirun -np 4 $QECONVERSE < Cl2x.in > Cl2x.out
mpirun -np 4 $QECONVERSE < Cl2y.in > Cl2y.out
mpirun -np 4 $QECONVERSE < Cl2z.in > Cl2z.out

python3 ../check_nmr.py \
    --files Na1x.out Na1y.out Na1z.out Cl2x.out Cl2y.out Cl2z.out \
    --refdir reference
