#ifndef IC_H_
#define IC_H_

#include "ctx.h"

PetscErrorCode set_ic_aug( Ctx *, Vec );
PetscErrorCode set_ic_aug_from_file( const char *, Ctx *, Vec );
PetscErrorCode set_ic_dSdr( Ctx *, Vec );

#endif
