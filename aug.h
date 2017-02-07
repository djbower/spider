#ifndef AUG_H_
#define AUG_H_

PetscErrorCode CreateAug(Vec in, Vec *out_aug);
PetscErrorCode CreateUnAug(Vec in_aug, Vec *out);
PetscErrorCode ToAug(Vec in,Vec out_aug);
PetscErrorCode FromAug(Vec in_aug,Vec out);

#endif
