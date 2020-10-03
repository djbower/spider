#if !defined(EOS_OUTPUT_H_)
#define EOS_OUTPUT_H_

/* All for outputting EOS-related data should be here */

#include "ctx.h"
#include "eos.h"

PetscErrorCode JSON_add_phase_boundary(const Ctx*,const EOS,const char*,cJSON*);

#endif
