#ifndef EOS_H_
#define EOS_H_

#include "ctx.h"

PetscErrorCode set_eos( Parameters * );

/* TODO: below probably moves elsewhere eventually (becomes static?) */
PetscErrorCode set_rtpress_struct( PetscScalar, PetscScalar, Ctx * );

/* for function testing */
PetscScalar get_rtpress_pressure_test( Ctx * );
PetscScalar get_rtpress_entropy_test( Ctx * );

/* to free memory */
void free_interp1d( Interp1d * );
void free_interp2d( Interp2d * );

PetscScalar get_val1d( Interp1d const *, PetscScalar );
PetscScalar get_val2d( Interp2d const *, PetscScalar, PetscScalar );

#endif
