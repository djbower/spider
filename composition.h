#ifndef COMPOSITION_H_
#define COMPOSITION_H_

#include "ctx.h"

PetscErrorCode set_rheological_front( Ctx * );
PetscErrorCode set_magma_ocean_crystal_fraction( Ctx *E );
PetscErrorCode set_magma_ocean_bridgmanite_fraction( Ctx *E );

#endif
