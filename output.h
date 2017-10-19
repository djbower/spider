#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "ctx.h"

PetscErrorCode scale_vectors_and_output( Vec *, PetscScalar *, PetscInt, PetscViewer );
PetscErrorCode atmosphere_structs_to_vec( Vec, Ctx *, Vec );
PetscErrorCode constants_struct_to_vec( Constants const *, Vec );
PetscErrorCode SetScalingsForOutput( Ctx * );
#endif
