#ifndef AUG_H_
#define AUG_H_

#include "petsc.h"

PetscErrorCode CreateAug(Vec in, Vec *out_aug);
PetscErrorCode CreateUnAug(Vec in_aug, Vec *out);
PetscErrorCode ToAug(Vec in,Vec out_aug);
PetscErrorCode FromAug(Vec in_aug,Vec out);

#endif
