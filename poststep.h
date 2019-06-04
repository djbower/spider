#if !defined(POSTSTEP_H_)
#define POSTSTEP_H_
#include "ctx.h"
#include <petscts.h>

typedef struct {
  PetscReal ref_concentration;
} PostStepData;

PetscErrorCode PostStepDataInitialize(Ctx *,Vec);
PetscErrorCode PostStep(TS);

#endif
