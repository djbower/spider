#if !defined(EOS_ADAMSWILLIAMSON_H_)
#define EOS_ADAMSWILLIAMSON_H_

#include "eos.h"

typedef struct {
    PetscScalar radius;
    PetscScalar radius_core;
    PetscScalar gravity;
    PetscScalar density_surface;
    PetscScalar density_average; // computed
    PetscScalar beta;
} data_EOSAdamsWilliamson;

#endif

PetscErrorCode EOSAdamsWilliamson_ObjectiveFunctionRadius( SNES, Vec, Vec, void* );
PetscErrorCode EOSAdamsWilliamson_JacobianRadius( SNES, Vec, Mat, Mat, void* );
PetscErrorCode EOSAdamsWilliamsonMassCoordinateSpatialDerivative( EOS eos, PetscScalar, PetscScalar, PetscScalar * );
