#ifndef UTIL_H_
#define UTIL_H_

#include "ctx.h"

PetscScalar average( PetscScalar, PetscScalar );
PetscScalar combine_matprop( PetscScalar, PetscScalar, PetscScalar );
PetscErrorCode d_dr( Ctx *, Vec, Vec );

#endif
