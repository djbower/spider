#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "ctx.h"

PetscErrorCode add_vector_to_binary_output( Vec, PetscViewer );
PetscErrorCode atmosphere_structs_to_vec( Atmosphere const *,AtmosphereParameters const *, Vec );

#endif
