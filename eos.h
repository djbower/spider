#ifndef EOS_H_
#define EOS_H_

#include "ctx.h"
#include "parameters.h"

PetscErrorCode EosParametersCreate( EosParameters * );
PetscErrorCode EosParametersDestroy( EosParameters * );
PetscErrorCode EosParametersSetFromOptions( EosParameters, const FundamentalConstants, const ScalingConstants );

/* TODO: below probably moves elsewhere eventually (becomes static?) */
PetscErrorCode set_rtpress_struct( PetscScalar, PetscScalar, Ctx * );

#if 0
/* TODO: need to refresh with new eos evaluation structs */
/* for function testing */
PetscScalar get_rtpress_pressure_test( Ctx * );
PetscScalar get_rtpress_entropy_test( Ctx * );
#endif

PetscScalar get_val1d( Interp1d const, PetscScalar );
PetscScalar get_val2d( Interp2d const, PetscScalar, PetscScalar );

#endif
