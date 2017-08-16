#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "ctx.h"

PetscScalar get_initial_xCO2( Atmosphere * );
PetscScalar get_initial_xH2O( Atmosphere * );
PetscErrorCode set_emissivity( Atmosphere * );
PetscErrorCode set_dxdt( Ctx * );

#endif
