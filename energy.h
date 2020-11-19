#ifndef ENERGY_H_
#define ENERGY_H_

#include "ctx.h"

PetscScalar GetConvectiveHeatFlux( Ctx *, PetscInt * );
PetscScalar GetConductiveHeatFlux( Ctx *, PetscInt * );
PetscScalar GetMixingHeatFlux( Ctx *, PetscInt *);
PetscScalar GetGravitationalHeatFlux( Ctx *, PetscInt * );
PetscErrorCode set_Etot( Ctx * );
PetscErrorCode set_Htot( Ctx *, PetscReal t );
PetscErrorCode set_interior_structure_from_solution( Ctx *, PetscScalar, Vec );
PetscScalar tsurf_param( PetscScalar, const AtmosphereParameters * );

#endif
