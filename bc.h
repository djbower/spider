#ifndef BC_H_
#define BC_H_

#include "ctx.h"

PetscErrorCode set_core_mantle_flux( Ctx * );
PetscErrorCode set_surface_flux( Ctx * );

PetscScalar get_initial_volatile( AtmosphereParameters const *, VolatileParameters const *, PetscScalar );
PetscErrorCode set_atmosphere_volatile_content( Ctx *, PetscScalar, PetscScalar );
PetscScalar get_dx0dt( Atmosphere *, AtmosphereParameters const *, PetscScalar, PetscScalar );
PetscScalar get_dx1dt( Atmosphere *, AtmosphereParameters const *, PetscScalar, PetscScalar );

#endif
