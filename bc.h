#ifndef AUG_H_
#define AUG_H_

#include "ctx.h"

PetscErrorCode get_core_cooling( Ctx * );
PetscScalar radiative_flux_with_dT( PetscScalar );

#endif
