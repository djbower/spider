#ifndef COMPOSITION_H_
#define COMPOSITION_H_

#include "petsc.h"
//#include "parameters.h"

typedef struct Composition_ {
    PetscScalar Brg_fraction; // molar fraction
    /* molar fraction of residual is res_fraction = 1 - brg_fraction */
    PetscScalar BSE_Brg_mass_ratio;
    /* above is kg/mol / kg/mol */
} Composition;

PetscScalar get_BSE_Brg_mass_ratio( PetscScalar, PetscScalar );

/* FIXME: below */
#if 0
PetscErrorCode set_composition( Ctx * );
PetscErrorCode initialise_composition( Ctx * );
#endif

#endif
