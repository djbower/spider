#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "ctx.h"

PetscErrorCode set_initial_xCO2( Atmosphere * );
PetscErrorCode set_initial_xH2O( Atmosphere * );
PetscErrorCode set_emissivity_abe_matsui( Atmosphere * );
PetscErrorCode set_dxdt( Ctx * );

#endif
