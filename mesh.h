#include "ctx.h"

PetscErrorCode set_mesh( Ctx * );
static PetscErrorCode spherical_area( DM, Vec, Vec );
static PetscErrorCode spherical_volume( Ctx *, Vec, Vec );
static PetscErrorCode mixing_length( DM, Vec, Vec );
static PetscErrorCode aw_pressure( DM, Vec, Vec );
static PetscErrorCode aw_pressure_gradient( DM, Vec, Vec );
