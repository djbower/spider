#ifndef BC_H_
#define BC_H_

#include "ctx.h"

PetscErrorCode set_core_mantle_flux( Ctx * );
PetscErrorCode set_surface_flux( Ctx * );

// atmosphere
PetscScalar get_initial_volatile( Ctx const *, VolatileParameters const * );
PetscErrorCode set_atmosphere_volatile_content( Ctx *, PetscScalar, PetscScalar );
PetscScalar get_dx0dt( Ctx const *, PetscScalar );
PetscScalar get_dx1dt( Ctx const *, PetscScalar );

#endif
