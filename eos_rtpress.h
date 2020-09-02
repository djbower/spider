#if !defined(EOS_RTPRESS_H_)
#define EOS_RTPRESS_H_

#include "petsc.h"
#include "eos.h"

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
    PetscScalar const * AVOGADRO_ptr;
} data_EOSRTpress;
typedef data_EOSRTpress* EOSRTpress;


// TODO not clear if these ultimately need to be in this header, or can just be static functions in eos_rtpress.c
PetscErrorCode RTpressParametersCreate( RTpressParameters * );
PetscErrorCode RTpressParametersCreateAndSet( RTpressParameters *, const FundamentalConstants );
PetscErrorCode RTpressParametersDestroy( RTpressParameters * );
PetscErrorCode SetEosEvalFromRTpress( const RTpressParameters, PetscScalar, PetscScalar, EosEval * );

#endif
