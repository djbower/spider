#if !defined(CONSTANTS_H_)
#define CONSTANTS_H_

#include "petsc.h"

/*
 ******************************************************************************
 * Dimensional constants
 ******************************************************************************
 */

/* constants to scale the physical problem, largely chosen based on numerical
   considerations */
typedef struct
{
    /* primary */
    PetscScalar RADIUS;
    PetscScalar TEMP;
    PetscScalar ENTROPY; /* note: specific entropy */
    PetscScalar DENSITY;
    PetscScalar VOLATILE;
    /* derived from primary */
    PetscScalar AREA;
    PetscScalar VOLUME;
    PetscScalar MASS;
    PetscScalar TIME;
    PetscScalar TIMEYRS;
    PetscScalar SENERGY; /* specific energy */
    PetscScalar ENERGY;
    PetscScalar PRESSURE;
    PetscScalar POWER;
    PetscScalar FLUX;
    PetscScalar DPDR;
    PetscScalar GRAVITY;
    PetscScalar KAPPA;
    PetscScalar DSDP;
    PetscScalar DSDR;
    PetscScalar DTDP;
    PetscScalar DTDR;
    PetscScalar GSUPER;
    PetscScalar VISC;
    PetscScalar LOG10VISC;
    PetscScalar COND;
    PetscScalar SIGMA;
    PetscScalar RHS;
    PetscScalar HEATGEN;
} data_ScalingConstants;
typedef data_ScalingConstants *ScalingConstants;

PetscErrorCode ScalingConstantsCreate(ScalingConstants *);
PetscErrorCode ScalingConstantsDestroy(ScalingConstants *);

/* fundamental constants */
typedef struct
{
    PetscScalar AVOGADRO;
    PetscScalar BOLTZMANN;
    PetscScalar GAS;
    PetscScalar GRAVITATIONAL;
    PetscScalar STEFAN_BOLTZMANN;
    PetscScalar OCEAN_MOLES;
} data_FundamentalConstants;
typedef data_FundamentalConstants *FundamentalConstants;

PetscErrorCode FundamentalConstantsCreate(FundamentalConstants *);
PetscErrorCode FundamentalConstantsDestroy(FundamentalConstants *);
PetscErrorCode FundamentalConstantsSet(FundamentalConstants, ScalingConstants const);

#endif
