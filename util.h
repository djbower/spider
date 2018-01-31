#ifndef UTIL_H_
#define UTIL_H_

#include "ctx.h"

PetscScalar combine_matprop( PetscScalar, PetscScalar, PetscScalar );
PetscScalar tanh_weight( PetscScalar, PetscScalar, PetscScalar );
PetscErrorCode set_d_dr( Ctx * );
PetscErrorCode set_entropy( Ctx *, PetscScalar );
PetscErrorCode set_dSdr_b_from_S_s( Ctx * );

#endif
