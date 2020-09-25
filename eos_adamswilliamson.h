#if !defined(EOS_ADAMSWILLIAMSON_H_)
#define EOS_ADAMSWILLIAMSON_H_

#include "eos.h"

typedef struct {
    PetscScalar radius;
    PetscScalar gravity;
    PetscScalar rhos;
    PetscScalar beta;
} data_EOSAdamsWilliamson;

PetscErrorCode EOSAdamsWilliamson_GetRadiusFromPressure( EOS, PetscScalar, PetscScalar * );
PetscErrorCode EOSAdamsWilliamson_GetMassWithinPressure( EOS, PetscScalar, PetscScalar, PetscScalar *);

#endif
