#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "ctx.h"

PetscErrorCode add_vector_to_binary_output( Vec, PetscViewer );
PetscErrorCode atmosphere_struct_to_vec( Atmosphere *, Vec );

#endif