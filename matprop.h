#ifndef MATPROP_H_
#define MATPROP_H_

#include "ctx.h"

PetscErrorCode set_melt_fraction_staggered( Ctx * );
PetscErrorCode set_capacitance_staggered( Ctx * );
PetscErrorCode set_matprop_basic( Ctx * );

#endif
