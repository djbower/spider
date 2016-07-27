#include "ctx.h"

PetscErrorCode set_core_cooling( Ctx * );
static PetscScalar utbl_temp_drop( PetscScalar );
PetscScalar radiative_flux_with_dT( PetscScalar );
static PetscScalar radiative_flux( PetscScalar );
