#if !defined(POSTSTEP_H_)
#define POSTSTEP_H_
#include "ctx.h"
#include "parameters.h"
#include <petscts.h>

typedef struct {
  PetscReal volatile_abundances[SPIDER_MAX_VOLATILE_SPECIES];
  PetscReal tsurf;
} PostStepData;

PetscErrorCode PostStepDataInitialize(Ctx *,Vec);
PetscErrorCode PostStep(TS);

#endif
