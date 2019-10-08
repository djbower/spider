#ifndef BC_H_
#define BC_H_

#include "ctx.h"

PetscErrorCode set_core_mantle_flux( Ctx * );
PetscErrorCode set_surface_flux( Ctx * );
PetscErrorCode solve_dxdts( Ctx * );
PetscScalar tsurf_param( PetscScalar, const AtmosphereParameters * );

#endif
