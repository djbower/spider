#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "parameters.h"

typedef struct Atmosphere_ {
    // calculated quantities (18)
    PetscScalar M0; // total mass of mantle from EOS (kg)
    PetscScalar Mliq; // mass of liquid (kg)
    PetscScalar Msol; // mass of solid (kg)
    PetscScalar dMliqdt; // dMliq/dt (kg/yr)
    PetscScalar tau; // aggregate optical depth (dimensionless)
    PetscScalar x0; // CO2 content (wt %)
    PetscScalar dx0dt; // dx0/dt (wt % / yr)
    PetscScalar p0; // CO2 partial pressure (Pa)
    PetscScalar dp0dx; // dp0/dx (Pa/mass fraction)
    PetscScalar m0; // CO2 mass in atmosphere (kg)
    PetscScalar tau0; // CO2 optical depth (dimensionless)
    PetscScalar x1; // H2O content (wt %)
    PetscScalar dx1dt; // dx1/dt (wt % / yr)
    PetscScalar p1; // H2O partial pressure (Pa)
    PetscScalar dp1dx; // dp1dx (Pa / mass fraction)
    PetscScalar m1; // H2O mass in atmosphere (kg)
    PetscScalar tau1; // H20 optical depth (dimensionless)
    PetscScalar emissivity; // variable emissivity (see also EMISSIVITY0 in AtmosphereParameters)
} Atmosphere;

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

PetscErrorCode set_initial_xCO2( Atmosphere *, AtmosphereParameters const * );
PetscErrorCode set_initial_xH2O( Atmosphere *, AtmosphereParameters const * );
PetscErrorCode set_emissivity_abe_matsui( Atmosphere *, AtmosphereParameters const * );
PetscErrorCode set_dx0dt( Atmosphere *, AtmosphereParameters const *);
PetscErrorCode set_dx1dt( Atmosphere *, AtmosphereParameters const *);

#endif
