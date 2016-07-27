#include "ctx.h"

PetscErrorCode set_twophase( Ctx * );
static PetscErrorCode set_liquidus( Ctx * );
static PetscErrorCode set_solidus( Ctx * );
static PetscErrorCode set_fusion( Ctx * );
static PetscErrorCode set_fusion_curve( Ctx * );
static PetscErrorCode set_mixed_phase( Ctx * );
