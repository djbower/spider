#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <petsc.h>
#include "lookup.h"

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
    PetscScalar VOLSCALE;
} Constants;

/* atmosphere */
typedef struct Atmosphere_ {
    // calculated quantities (14)
    PetscScalar Mliq; // mass of liquid (kg)
    PetscScalar Msol; // mass of solid (kg)
    PetscScalar dMliqdt; // dMliq/dt (kg/yr)
    PetscScalar tsurf; // surface temperature
    PetscScalar tau; // aggregate optical depth (dimensionless)
    PetscScalar p0; // CO2 partial pressure (Pa)
    PetscScalar dp0dx; // dp0/dx (Pa/mass fraction)
    PetscScalar m0; // CO2 mass in atmosphere (kg)
    PetscScalar tau0; // CO2 optical depth (dimensionless)
    PetscScalar p1; // H2O partial pressure (Pa)
    PetscScalar dp1dx; // dp1dx (Pa / mass fraction)
    PetscScalar m1; // H2O mass in atmosphere (kg)
    PetscScalar tau1; // H20 optical depth (dimensionless)
    PetscScalar emissivity; // variable emissivity (see also EMISSIVITY0 in AtmosphereParameters)
} Atmosphere;

typedef struct VolatileParameters_ {
    PetscScalar initial;
    PetscScalar kdist;
    PetscScalar kabs;
    PetscScalar henry;
    PetscScalar henry_pow;
} VolatileParameters;

/* for storing atmosphere outputs for eventual writing to Petsc
   binary file */
typedef enum {MO_ATMOSPHERE_TYPE_GREY_BODY=1,MO_ATMOSPHERE_TYPE_ZAHNLE,MO_ATMOSPHERE_TYPE_VOLATILES} MagmaOceanAtmosphereType;
typedef struct AtmosphereParameters_ {
    // input parameters
    MagmaOceanAtmosphereType MODEL;
    PetscBool HYBRID;
    // below are standard, also used for grey-body atmosphere
    PetscScalar emissivity0;
    PetscScalar sigma;
    PetscScalar teqm;
    PetscBool   PARAM_UTBL;
    PetscScalar param_utbl_const;
    // for volatile ODE
    PetscBool SOLVE_FOR_VOLATILES;
    PetscScalar P0; 
    VolatileParameters H2O_volatile_parameters;
    VolatileParameters CO2_volatile_parameters;
} AtmosphereParameters;

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
    char        liquidusFilename[PETSC_MAX_PATH_LEN];
    char        solidusFilename[PETSC_MAX_PATH_LEN];
    char        alphaSolFilename[PETSC_MAX_PATH_LEN];
    char        alphaMelFilename[PETSC_MAX_PATH_LEN];
    char        cpSolFilename[PETSC_MAX_PATH_LEN];
    char        cpMelFilename[PETSC_MAX_PATH_LEN];
    char        dtdpsSolFilename[PETSC_MAX_PATH_LEN];
    char        dtdpsMelFilename[PETSC_MAX_PATH_LEN];
    char        rhoSolFilename[PETSC_MAX_PATH_LEN];
    char        rhoMelFilename[PETSC_MAX_PATH_LEN];
    char        tempSolFilename[PETSC_MAX_PATH_LEN];
    char        tempMelFilename[PETSC_MAX_PATH_LEN];
    Lookup      melt_prop;
    Lookup      solid_prop;

    // Additional Atmosphere Parameters
    AtmosphereParameters atmosphere_parameters;

    // Scaling factors / dimensionalising constants
    Constants constants;
} Parameters;

PetscErrorCode InitializeParametersAndSetFromOptions(Parameters *parameters);
PetscErrorCode PrintParameters(Parameters const *parameters);
PetscErrorCode SetLookups( Parameters * );

#endif
