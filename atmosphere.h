#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "ctx.h"

PetscErrorCode set_initial_carbon( Ctx * );
PetscErrorCode set_initial_water( Ctx * );
PetscErrorCode atmosphere_test( Ctx * );

#endif
