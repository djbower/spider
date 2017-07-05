#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "ctx.h"

PetscErrorCode set_initial_carbon( Ctx * );
PetscErrorCode set_initial_water( Ctx * );
PetscScalar get_emissivity( Ctx * );

#endif
