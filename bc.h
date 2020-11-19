#ifndef BC_H_
#define BC_H_

#include "ctx.h"

PetscErrorCode SetCoreMantleFluxBC( Ctx * );
PetscErrorCode set_surface_flux( Ctx * );
PetscErrorCode solve_dpdts( Ctx * );

#endif
