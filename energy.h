#ifndef ENERGY_H_
#define ENERGY_H_

#include "ctx.h"

PetscScalar GetConvectiveHeatFlux( Ctx *, PetscInt * );
PetscScalar GetConductiveHeatFlux( Ctx *, PetscInt * );
PetscScalar GetMixingHeatFlux( Ctx *, PetscInt *);
PetscScalar GetGravitationalHeatFlux( Ctx *, PetscInt * );
PetscErrorCode set_Etot( Ctx * );
PetscErrorCode set_Htot( Ctx *, PetscReal t );
PetscErrorCode set_current_state_from_solution( Ctx *, PetscScalar, Vec );
PetscErrorCode set_interior_atmosphere_interface_from_surface_entropy( Ctx * );
PetscErrorCode solve_for_surface_radiation_balance( Ctx *, PetscReal );

#endif
