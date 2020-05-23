#ifndef EOS_H_
#define EOS_H_

#include "ctx.h"
#include "parameters.h"

typedef struct RTpressEval_ {
  PetscScalar P; /* pressure */
  PetscScalar S; /* entropy */
  RTpressParameters rtp;
} RTpressEval;

PetscErrorCode EosParametersCreate( EosParameters * );
PetscErrorCode EosParametersDestroy( EosParameters * );
PetscErrorCode EosParametersSetFromOptions( EosParameters, const FundamentalConstants, const ScalingConstants );

PetscErrorCode SetEosEval( const EosParameters, PetscScalar, PetscScalar, EosEval * );
PetscErrorCode SetEosCompositeEval( const EosComposite, PetscScalar, PetscScalar, EosEval * );

PetscErrorCode EosCompositeCreateTwoPhase( EosComposite *, const EosParameters[], PetscInt );
PetscErrorCode EosCompositeDestroy( EosComposite * );

/* TODO: below probably moves elsewhere eventually (becomes static?) */
//PetscErrorCode set_rtpress_struct( PetscScalar, PetscScalar, Ctx * );

#if 0
/* TODO: need to refresh with new eos evaluation structs */
/* for function testing */
PetscScalar get_rtpress_pressure_test( Ctx * );
PetscScalar get_rtpress_entropy_test( Ctx * );
#endif

PetscErrorCode SetInterp1dValue( Interp1d const, PetscScalar, PetscScalar *, PetscScalar * );
PetscScalar GetInterp2dValue( Interp2d const, PetscScalar, PetscScalar );

#endif
