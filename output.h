#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "ctx.h"

PetscErrorCode add_vector_to_binary_output( Vec, PetscViewer );
PetscErrorCode atmosphere_structs_to_vec( Ctx *, Vec );
PetscErrorCode constants_struct_to_vec( Constants const *, Vec );
#endif
