#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "ctx.h"

PetscScalar get_initial_xCO2( Ctx * );
PetscScalar get_emissivity( Ctx *, PetscScalar, PetscScalar );
PetscScalar get_dX0dt( Ctx *, PetscScalar, Vec );

#endif
