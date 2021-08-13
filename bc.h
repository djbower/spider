#ifndef BC_H_
#define BC_H_

#include "ctx.h"

PetscErrorCode set_surface_entropy_gradient_update( Ctx *, Vec );
PetscErrorCode set_cmb_entropy_gradient_update( Ctx *, Vec );
PetscErrorCode solve_dpdts( Ctx * );
PetscErrorCode get_dpdts_from_lookup( Ctx * );
PetscErrorCode set_surface_entropy_from_surface_gradient( Ctx * );
PetscErrorCode set_cmb_entropy_constant( Ctx * );
PetscErrorCode set_surface_entropy_constant( Ctx * );
PetscErrorCode set_cmb_entropy_from_cmb_gradient( Ctx * );
PetscErrorCode set_surface_flux_from_atmosphere( Ctx * );
PetscScalar get_tsurf_using_parameterised_boundary_layer( PetscScalar, const AtmosphereParameters );
PetscScalar get_dtsurf_using_parameterised_boundary_layer( PetscScalar, const AtmosphereParameters );

#endif
