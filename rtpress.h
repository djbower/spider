#ifndef RTPRESS_H_
#define RTPRESS_H_

#include "ctx.h"

PetscErrorCode set_rtpress_parameters( Eos * );
PetscScalar get_rtpress_temperature( PetscScalar, PetscScalar, Ctx * );

#endif
