#ifndef RHEOLOGICALFRONT_H_
#define RHEOLOGICALFRONT_H_

#include <petscdm.h>
#include "parameters.h"
#include "cJSON.h"

typedef struct RheologicalFrontMantleProperties_ {
    PetscScalar phi;
    PetscScalar depth;
    PetscScalar pressure;
    PetscScalar temperature;
} RheologicalFrontMantleProperties;

typedef struct RheologicalFront_ {
    PetscInt mesh_index;
    PetscScalar depth;
    PetscScalar pressure;
    PetscScalar temperature;
    PetscScalar phi_global;
    /* next are different ways of computing the average mantle
       properties above and below the rheological front
         - middle uses the mid-point from the rheological front
           and either the CMB or surface (as appropriate)
         - mass_avg uses the mass-averaged quantity, again, for
           above and below the rheological front */
    RheologicalFrontMantleProperties above_middle;
    RheologicalFrontMantleProperties above_mass_avg;
    RheologicalFrontMantleProperties below_middle;
    RheologicalFrontMantleProperties below_mass_avg;
} RheologicalFront;

PetscErrorCode JSON_add_rheological_front( DM, ScalingConstants, RheologicalFront *, const char *, cJSON * );
PetscInt get_crossover_index( DM, Vec, PetscScalar, PetscInt );
#endif
