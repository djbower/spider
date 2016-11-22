#ifndef UTIL_H_
#define UTIL_H_

#include "ctx.h"

PetscScalar combine_matprop( PetscScalar, PetscScalar, PetscScalar );
PetscScalar tanh_weight( PetscScalar, PetscScalar, PetscScalar );
PetscErrorCode set_d_dr( Ctx * );
PetscErrorCode set_d_dr_linear( Ctx * );
PetscErrorCode set_d_dr_quadratic( Ctx * );

#endif
