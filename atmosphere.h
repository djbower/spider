#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "ctx.h"

PetscScalar get_initial_xCO2( Ctx * );
PetscErrorCode set_emissivity( Ctx * );
PetscErrorCode set_dx0dt( Ctx * );

#endif
