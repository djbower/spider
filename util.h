#ifndef UTIL_H_
#define UTIL_H_

#include "ctx.h"

PetscScalar combine_matprop( PetscScalar, PetscScalar, PetscScalar );
PetscScalar tanh_weight( PetscScalar, PetscScalar, PetscScalar );
PetscErrorCode set_d_dr( Ctx * );
PetscErrorCode set_entropy_from_solution( Ctx *, Vec );
PetscErrorCode set_solution_from_entropy( Ctx *, Vec );

#endif
