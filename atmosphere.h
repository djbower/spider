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
    PetscScalar volscale;
    PetscScalar P0;
    VolatileParameters H2O_volatile_parameters;
    VolatileParameters CO2_volatile_parameters;
    /* although RADIUS and GRAVITY are duplicated, they may have
       different non-dimensional values depending on the volatile
       non-dimensionalisation scheme, which will be different to
       the dS/dr scheme to ensure comparable magnitude residuals */
    PetscScalar RADIUS; // duplicate
    PetscScalar GRAVITY; // duplicate
} AtmosphereParameters;

PetscScalar get_initial_volatile( Atmosphere *, AtmosphereParameters const *, VolatileParameters const * );
PetscErrorCode set_emissivity_abe_matsui( Atmosphere *, AtmosphereParameters const * );
PetscErrorCode set_dx0dt( Atmosphere *, AtmosphereParameters const *);
PetscErrorCode set_dx1dt( Atmosphere *, AtmosphereParameters const *);

PetscScalar tsurf_param( PetscScalar, AtmosphereParameters const * );
PetscScalar grey_body( PetscScalar, Atmosphere *, AtmosphereParameters const * );
PetscScalar steam_atmosphere_zahnle_1988( PetscScalar, PetscScalar, PetscScalar );

#endif
