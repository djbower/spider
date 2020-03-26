#ifndef RTPRESS_H_
#define RTPRESS_H_

#include "ctx.h"

/* rtpress initialisation */
PetscErrorCode set_rtpress_parameters( Eos * );
PetscErrorCode set_rtpress_struct( PetscScalar, PetscScalar, Ctx * );

#endif
