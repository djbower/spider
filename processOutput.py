#!/usr/bin/env python

# This is a test/demonstration of reading PETSc binary files

# Import PETSc's python reader 
import os
import sys
PETSC_DIR=os.getenv('PETSC_DIR')
if not PETSC_DIR :
    print('You must define PETSC_DIR in your environment')
    sys.exit(1)
sys.path.append(os.path.join(PETSC_DIR,'bin'))
import PetscBinaryIO as pio

io = pio.PetscBinaryIO('double'); # We expect *double* binary output from our quadruple-precision build (achieved with a flag in PETSc)

# Pick a set of steps to import
stepsToImport = [0,3,7]

# How many vectors we expect in each output file
nVecs = 3

# Attempt to read from step files
for step in stepsToImport :
    print('#### STEP ' + str(step) + ' ########')
    fh = open(os.path.join('output','petscbin.') + str(step))
    for i in range(1,nVecs+1) :
        objecttype = io.readObjectType(fh)
        if objecttype == 'Vec':
            v = io.readVec(fh)
            print('>> Vector ' + str(i) + ' of ' + str(nVecs) + ':')
            print(v)
        else :
            print('Unexpected object type')
