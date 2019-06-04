#if !defined(POSTSTEP_H_)
#define POSTSTEP_H_
#include "ctx.h"
#include <petscts.h>

typedef struct {
  PetscReal x_CO2; // CO2 concentration in melt
  PetscReal x_H2O; // H2O concentration in melt
} PostStepData;

PetscErrorCode PostStepDataInitialize(Ctx *,Vec);
PetscErrorCode PostStep(TS);

#endif
