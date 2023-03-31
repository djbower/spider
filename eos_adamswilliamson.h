#if !defined(EOS_ADAMSWILLIAMSON_H_)
#define EOS_ADAMSWILLIAMSON_H_

#include "eos.h"

typedef struct
{
    PetscScalar radius;
    PetscScalar radius_core;
    PetscScalar gravity;
    PetscScalar density_surface;
    PetscScalar density_average; // computed
    PetscScalar beta;
} data_EOSAdamsWilliamson;

#endif

PetscErrorCode EOSAdamsWilliamsonGetMassWithinShell(const EOS, PetscScalar, PetscScalar, PetscScalar *);
PetscErrorCode EOSAdamsWilliamsonGetPressureFromRadius(const EOS, PetscScalar, PetscScalar *);
PetscErrorCode EOSAdamsWilliamsonGetPressureGradientFromRadius(const EOS, PetscScalar, PetscScalar *);
PetscErrorCode EOSAdamsWilliamson_ObjectiveFunctionRadius(SNES, Vec, Vec, void *);
PetscErrorCode EOSAdamsWilliamson_JacobianRadius(SNES, Vec, Mat, Mat, void *);
PetscErrorCode EOSAdamsWilliamsonMassCoordinateSpatialDerivative(const EOS, PetscScalar, PetscScalar, PetscScalar *);
