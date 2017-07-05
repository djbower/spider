#ifndef BC_H_
#define BC_H_

#include "ctx.h"

PetscErrorCode set_core_mantle_flux( Ctx * );
PetscErrorCode set_surface_flux( Ctx *, PetscReal, PetscScalar );

#endif
