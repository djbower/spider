#ifndef TWOPHASE_H_
#define TWOPHASE_H_

#include "ctx.h"

PetscErrorCode set_twophase( Ctx * );
PetscErrorCode set_gphi_smooth( Ctx * );
PetscErrorCode set_Mliq( Ctx * );
PetscErrorCode set_Msol( Ctx * );
PetscErrorCode set_dMliqdt( Ctx * );

PetscErrorCode set_rheological_front_phi( Ctx * );
PetscErrorCode add_rheological_front_to_cJSON( Ctx const *, RheologicalFront const *, cJSON * );

#endif
