#ifndef BC_H_
#define BC_H_

#include "ctx.h"

PetscErrorCode set_core_mantle_flux( Ctx * );
PetscErrorCode set_surface_flux( Ctx * );
// TODO: remove eventually
//PetscScalar get_dxdt( const Ctx *, const VolatileParameters *, const Volatile * );
//PetscScalar get_dx0dt( const Ctx *, PetscScalar );
//PetscScalar get_dx1dt( const Ctx *, PetscScalar );

#endif
