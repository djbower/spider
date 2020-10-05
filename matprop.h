#ifndef MATPROP_H_
#define MATPROP_H_

#include "ctx.h"

PetscErrorCode set_phase_fraction_staggered( Ctx * );
PetscErrorCode set_capacitance_staggered( Ctx * );
PetscErrorCode set_matprop_basic( Ctx * );
PetscErrorCode GetEddyDiffusivity( const EOSEvalData, const Parameters, PetscScalar, PetscScalar, PetscScalar, PetscScalar *, PetscScalar *, PetscScalar * );

#endif
