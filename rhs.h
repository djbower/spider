#ifndef RHS_H_
#define RHS_H_

#include "ctx.h"

PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec S_in,Vec rhs_s,void *ptr);

#endif
