#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "ctx.h"
#include "dimensionalisablefield.h"

PetscErrorCode scale_vectors_and_output( DimensionalisableField* , PetscInt, PetscViewer );
PetscErrorCode atmosphere_structs_to_vec( Vec, Ctx *, Vec );
PetscErrorCode constants_struct_to_vec( Constants const *, Vec );

#endif
