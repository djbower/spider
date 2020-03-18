#if !defined(POSTSTEP_H_)
#define POSTSTEP_H_
#include "ctx.h"
#include "parameters.h"
#include <petscts.h>

typedef struct {
  PetscReal volatile_partialp[SPIDER_MAX_VOLATILE_SPECIES];
  PetscReal tsurf;
} PostStepData;

PetscErrorCode PostStepDataInitialize(Ctx *);
PetscErrorCode PostStep(TS);

#endif
