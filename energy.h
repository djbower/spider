#ifndef ENERGY_H_
#define ENERGY_H_

#include "ctx.h"

PetscErrorCode set_Etot( Ctx * );
PetscErrorCode set_Htot( Ctx *, PetscReal t );

#endif
