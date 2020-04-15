#ifndef EOS_H_
#define EOS_H_

#include "ctx.h"
#include "parameters.h"

PetscErrorCode EosParametersCreate( EosParameters * );
PetscErrorCode EosParametersDestroy( EosParameters * );
PetscErrorCode EosParametersSetFromOptions( EosParameters, const FundamentalConstants, const ScalingConstants );

/* TODO: need equivalent function for RTpress */
PetscErrorCode SetEosEvalFromLookup( const Lookup, PetscScalar, PetscScalar, EosEval * );

/* TODO: below probably moves elsewhere eventually (becomes static?) */
//PetscErrorCode set_rtpress_struct( PetscScalar, PetscScalar, Ctx * );

#if 0
/* TODO: need to refresh with new eos evaluation structs */
/* for function testing */
PetscScalar get_rtpress_pressure_test( Ctx * );
PetscScalar get_rtpress_entropy_test( Ctx * );
#endif

PetscScalar GetInterp1dValue( Interp1d const, PetscScalar );
PetscScalar GetInterp2dValue( Interp2d const, PetscScalar, PetscScalar );

#endif
