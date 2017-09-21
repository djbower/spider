#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "parameters.h"

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

PetscScalar get_initial_volatile( AtmosphereParameters const *, VolatileParameters const *, PetscScalar );
PetscScalar get_emissivity_abe_matsui( Atmosphere *, AtmosphereParameters const * );
PetscScalar get_emissivity_from_flux( Atmosphere const *, AtmosphereParameters const *, PetscScalar );
PetscScalar get_dx0dt( Atmosphere *, AtmosphereParameters const *, PetscScalar, PetscScalar );
PetscScalar get_dx1dt( Atmosphere *, AtmosphereParameters const *, PetscScalar, PetscScalar );

PetscScalar tsurf_param( PetscScalar, AtmosphereParameters const * );
PetscScalar grey_body( Atmosphere const *, AtmosphereParameters const * );
PetscScalar steam_atmosphere_zahnle_1988( Atmosphere const *, PetscScalar, PetscScalar );
PetscErrorCode set_atmosphere_volatile_content( Atmosphere *, AtmosphereParameters const *, PetscScalar, PetscScalar );

#endif
