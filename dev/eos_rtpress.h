#if !defined(EOS_RTPRESS_H_)
#define EOS_RTPRESS_H_

#include "eos.h"
#include "petsc.h"

typedef struct {
    PetscScalar V0;
    PetscScalar T0;
    PetscScalar S0;
    PetscScalar K0;
    PetscScalar KP0;
    PetscScalar E0;
    PetscScalar gamma0;
    PetscScalar gammaP0;
    PetscScalar m;
    PetscScalar b0;
    PetscScalar b1;
    PetscScalar b2;
    PetscScalar b3;
    PetscScalar b4;
    PetscScalar mavg;
    PetscScalar PV_UNIT;
    PetscScalar KBOLTZ;
    PetscScalar bscale;
    PetscScalar AVOGADRO;
} data_EOSRTpress;

#endif
