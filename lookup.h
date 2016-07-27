#include "ctx.h"

PetscErrorCode set_lookups( Ctx * );
static PetscErrorCode set_interp2d( const char *, Interp2d * );
static PetscErrorCode set_interp1d( const char *, Interp1d *, PetscInt );
PetscScalar get_val1d( Interp1d *, PetscScalar );
PetscScalar get_val2d( Interp2d *, PetscScalar, PetscScalar );
PetscErrorCode free_interp1d( Interp1d * );
PetscErrorCode free_interp2d( Interp2d * );
PetscErrorCode free_memory_interp( Ctx * );
