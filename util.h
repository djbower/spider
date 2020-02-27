#ifndef UTIL_H_
#define UTIL_H_

#include "ctx.h"

PetscScalar* Make1DPetscScalarArray( PetscInt );
PetscScalar** Make2DPetscScalarArray( PetscInt, PetscInt );
PetscScalar combine_matprop( PetscScalar, PetscScalar, PetscScalar );
PetscScalar tanh_weight( PetscScalar, PetscScalar, PetscScalar );
PetscErrorCode set_d_dr( Ctx * );
PetscErrorCode set_entropy_from_solution( Ctx *, Vec );
PetscErrorCode set_solution_from_entropy_at_staggered_nodes( Ctx *, Vec );
PetscErrorCode set_volatile_abundances_from_solution( Ctx *, Vec );
PetscErrorCode set_solution_from_partial_pressures( Ctx *, Vec );
PetscErrorCode average_by_mass_staggered( Ctx *, Vec, Vec *, PetscScalar * );
PetscErrorCode invert_vec_mask( Vec * );
PetscErrorCode make_vec_mask( DM, PetscInt, Vec * );

#endif
