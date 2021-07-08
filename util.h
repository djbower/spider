#ifndef UTIL_H_
#define UTIL_H_

#include "ctx.h"

PetscErrorCode MakeRelativeToSourcePathAbsolute(char*);
PetscErrorCode Make2DPetscScalarArray( PetscInt, PetscInt, PetscScalar *** );
PetscScalar combine_matprop( PetscScalar, PetscScalar, PetscScalar );
PetscScalar tanh_weight( PetscScalar, PetscScalar, PetscScalar );
PetscScalar get_smoothing( PetscScalar, PetscScalar );
PetscErrorCode set_d_dxi( Ctx * );
PetscErrorCode set_entropy_from_solution( Ctx *, Vec );
PetscErrorCode set_entropy_reconstruction_from_ctx( Ctx *, PetscScalar );
PetscErrorCode set_solution_from_entropy( Ctx *, Vec );
PetscErrorCode set_partial_pressures_from_solution( Ctx *, Vec );
PetscErrorCode set_solution_from_partial_pressures( Ctx *, Vec );
PetscErrorCode average_by_mass_staggered( Ctx *, Vec, Vec *, PetscScalar * );
PetscErrorCode invert_vec_mask( Vec * );
PetscErrorCode make_vec_mask( DM, PetscInt, Vec * );
PetscErrorCode PetscScalarCheckPositive( PetscScalar, const char * );
PetscErrorCode PetscIntCheckPositive( PetscInt, const char * );
PetscErrorCode PetscOptionsGetPositiveScalar( const char *, PetscScalar *, PetscScalar, PetscBool * );
PetscErrorCode PetscOptionsGetPositiveInt( const char *, PetscInt *, PetscInt, PetscBool * );
#endif
