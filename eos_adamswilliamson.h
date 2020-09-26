#if !defined(EOS_ADAMSWILLIAMSON_H_)
#define EOS_ADAMSWILLIAMSON_H_

#include "eos.h"

typedef struct {
    PetscScalar radius;
    PetscScalar radius_core;
    PetscScalar gravity;
    PetscScalar rhos;
    PetscScalar beta;
} data_EOSAdamsWilliamson;

PetscErrorCode EOSAdamsWilliamson_GetMassCoordinateAverageRho( EOS, PetscScalar * );

#endif
