#ifndef LOOKUP_H_
#define LOOKUP_H_

#include "ctx.h"

PetscErrorCode set_lookups( Ctx * );
PetscScalar get_val1d( Interp1d *, PetscScalar );
PetscScalar get_val2d( Interp2d *, PetscScalar, PetscScalar );
PetscErrorCode free_interp1d( Interp1d * );
PetscErrorCode free_interp2d( Interp2d * );
PetscErrorCode free_memory_interp( Ctx * );

#endif
