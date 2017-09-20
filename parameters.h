#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <petsc.h>
#include "lookup.h"
#include "atmosphere.h"


/* dimensionalising constants */
typedef struct _Constants {
    // primary
    PetscScalar RADIUS;
    PetscScalar TEMP;
    PetscScalar ENTROPY;
    PetscScalar DENSITY;
    // derived from primary
    PetscScalar AREA;
    PetscScalar VOLUME;
    PetscScalar MASS;
    PetscScalar TIME;
    PetscScalar TIMEYRS;
    PetscScalar SENERGY;
    PetscScalar ENERGY;
    PetscScalar PRESSURE;
    PetscScalar POWER;
    PetscScalar FLUX;
    PetscScalar DPDR;
    PetscScalar GRAVITY;
    PetscScalar KAPPA;
    PetscScalar DTDP;
    PetscScalar DSDR;
    PetscScalar DTDR;
    PetscScalar GSUPER;
    PetscScalar VISC;
    PetscScalar LOG10VISC;
    PetscScalar COND;
    PetscScalar SIGMA;
    PetscScalar LHS;
    PetscScalar RHS;
} Constants;

typedef struct _Parameters {

    // Discretization parameters
    PetscInt    nstepsmacro,maxsteps;
    PetscReal   dtmacro;
    PetscReal   t0; /* Initial time */
    PetscInt    numpts_b,numpts_s;

    PetscBool monitor;

    //  "Standard" parameters
    // 19
    PetscScalar sinit;
    PetscScalar ic_dsdr;
    PetscScalar radius;
    PetscScalar coresize;
    PetscScalar rhos;
    PetscScalar beta;
    PetscScalar grain;
    PetscScalar gravity;
    PetscScalar phi_critical;
    PetscScalar phi_width;
    PetscScalar phi_skew;
    PetscScalar rho_core;
    PetscScalar cp_core;
    PetscScalar tfac_core_avg;
    PetscScalar swidth;
    PetscScalar log10visc_sol;
    PetscScalar cond_sol;
    PetscScalar log10visc_mel;
    PetscScalar cond_mel;

    // Lookup tables
    Lookup      melt_prop;
    Lookup      solid_prop;

    // Additional Atmosphere Parameters
    AtmosphereParameters atmosphere_parameters;

    // Scaling factors / dimensionalizing constants
    Constants constants;
} Parameters;

PetscErrorCode InitializeParametersAndSetFromOptions(Parameters *parameters);
PetscErrorCode PrintParameters(Parameters const *parameters);
PetscErrorCode SetLookups( Parameters * );

#endif
