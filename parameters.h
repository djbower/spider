#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <petsc.h>

/* for storing atmosphere outputs for eventual writing to Petsc
   binary file */
typedef enum {MO_ATMOSPHERE_TYPE_GREY_BODY=1,MO_ATMOSPHERE_TYPE_ZAHNLE,MO_ATMOSPHERE_TYPE_VOLATILES} MagmaOceanAtmosphereType;
typedef struct AtmosphereParameters_ {
    // input parameters (20)
    MagmaOceanAtmosphereType MODEL;
    PetscInt HYBRID;
    // below are standard, also used for grey-body atmosphere
    PetscScalar EMISSIVITY0;
    PetscScalar SIGMA;
    PetscScalar TEQM;
    PetscScalar CONSTBC;
    // for volatile ODE
    PetscScalar VOLSCALE;
    PetscScalar P0;
    // H2O (TODO: move to a 'Volatile' struct)
    PetscScalar H2O_INITIAL;
    PetscScalar H2O_KDIST;
    PetscScalar H2O_KABS;
    PetscScalar H2O_HENRY;
    PetscScalar H2O_HENRY_POW;
    // CO2 (TODO: move to a 'Volatile' struct)
    PetscScalar CO2_INITIAL;
    PetscScalar CO2_KDIST;
    PetscScalar CO2_KABS;
    PetscScalar CO2_HENRY;
    PetscScalar CO2_HENRY_POW;
    /* although RADIUS and GRAVITY are duplicated, they may have
       different non-dimensional values depending on the volatile
       non-dimensionalisation scheme, which will be different to
       the dS/dr scheme to ensure comparable magnitude residuals */
    PetscScalar RADIUS; // duplicate
    PetscScalar GRAVITY; // duplicate
} AtmosphereParameters;

/* dimensionalising constants */
typedef struct _Constants {
    PetscScalar RADIUS;
    PetscScalar TEMP;
    PetscScalar ENTROPY;
    PetscScalar DENSITY;
    PetscScalar AREA;
    PetscScalar AREAG; // with 4*pi geometry
    PetscScalar VOLUME;
    PetscScalar VOLUMEG; // with 4*pi geometry
    PetscScalar MASS;
    PetscScalar MASSG; // with 4*pi geometry
    PetscScalar TIME;
    PetscScalar TIMEYRS;
    PetscScalar SENERGY;
    PetscScalar ENERGY;
    PetscScalar ENERGYG; // with 4*pi geometry
    PetscScalar PRESSURE;
    PetscScalar POWER;
    PetscScalar POWERG; // with 4*pi geometry
    PetscScalar FLUX;
    PetscScalar DPDR;
    PetscScalar ALPHA;
    PetscScalar GRAVITY;
    PetscScalar KAPPA;
    PetscScalar DTDP;
    PetscScalar DSDR;
    PetscScalar DTDR;
    PetscScalar GSUPER;
    PetscScalar ETA;
    PetscScalar LOG10ETA;
    PetscScalar NU;
    PetscScalar COND;
    PetscScalar SIGMA;
    PetscScalar LHS;
    PetscScalar LHSG; // with 4*pi geometry
    PetscScalar RHS;
} Constants;

typedef struct _Parameters {
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

    // Additional Atmosphere Parameters
    AtmosphereParameters atmosphere_parameters;

    // Scaling factors / dimensionalizing constants
    Constants constants;
} Parameters;

PetscErrorCode InitializeParameters(Parameters *parameters);
PetscErrorCode SetParametersFromOptions(Parameters *parameters);
PetscErrorCode PrintParameters(Parameters const *parameters,FILE *file);

#endif
