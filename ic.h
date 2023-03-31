#ifndef IC_H_
#define IC_H_

#include "ctx.h"

PetscErrorCode set_initial_condition(Ctx *, Vec);
PetscErrorCode read_JSON_file_to_JSON_object(const char *, cJSON **);

#endif
