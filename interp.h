#ifndef INTERP_H_
#define INTERP_H_

#include <petsc.h>

typedef struct {
    PetscInt    NX;
    PetscScalar *xa;
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar *ya;
    PetscScalar ymin;
    PetscScalar ymax;
} data_Interp1d;
typedef data_Interp1d* Interp1d;

typedef struct {
    PetscInt    NX;
    PetscScalar *xa;
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar dx;
    PetscInt    NY;
    PetscScalar *ya;
    PetscScalar ymin;
    PetscScalar ymax;
    PetscScalar dy;
    PetscScalar **za;
} data_Interp2d;
typedef data_Interp2d* Interp2d;

PetscErrorCode Interp1dCreateAndSet( const char *, Interp1d *, PetscScalar, PetscScalar );
PetscErrorCode Interp1dDestroy( Interp1d * );
PetscErrorCode Interp2dCreateAndSet( const char *, Interp2d *, PetscScalar, PetscScalar, PetscScalar );
PetscErrorCode Interp2dDestroy( Interp2d * );
PetscErrorCode SetInterp1dValue( Interp1d const, PetscScalar, PetscScalar *, PetscScalar * );
PetscErrorCode SetInterp2dValue( Interp2d const, PetscScalar, PetscScalar, PetscScalar * );

#endif
