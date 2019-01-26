#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include <petsc.h>
#include "parameters.h"
#include "cJSON.h"

typedef struct Volatile_ {
    // FIXME: need to update this quantity from the RHS
    PetscScalar x; // ppm in liquid mantle
    PetscScalar p; // partial pressure (Pa)
    PetscScalar dpdx; // dp/dx (Pa/mass fraction)
    PetscScalar m; // mass in atmosphere (kg)
    PetscScalar tau; // optical_depth (non-dimensional)
} Volatile;

typedef struct Atmosphere_ {
    // calculated quantities (14)
    PetscScalar Mliq; // mass of liquid (kg)
    PetscScalar Msol; // mass of solid (kg)
    PetscScalar dMliqdt; // dMliq/dt (kg/yr)
    PetscScalar tsurf; // surface temperature
    PetscScalar tau; // aggregate optical depth (dimensionless)
    PetscScalar emissivity; // variable emissivity (see also EMISSIVITY0 in AtmosphereParameters)
    Volatile    CO2;
    Volatile    H2O;
    /* TODO: remove, now in Volatile struct */
    //PetscScalar p0; // CO2 partial pressure (Pa)
    //PetscScalar dp0dx; // dp0/dx (Pa/mass fraction)
    //PetscScalar m0; // CO2 mass in atmosphere (kg)
    //PetscScalar tau0; // CO2 optical depth (dimensionless)
    //PetscScalar p1; // H2O partial pressure (Pa)
    //PetscScalar dp1dx; // dp1dx (Pa / mass fraction)
    //PetscScalar m1; // H2O mass in atmosphere (kg)
    //PetscScalar tau1; // H20 optical depth (dimensionless)
} Atmosphere;

PetscScalar get_grey_body_flux( const Atmosphere *, const AtmosphereParameters * );
PetscScalar get_steam_atmosphere_zahnle_1988_flux( const Atmosphere *, const Constants *C );
PetscScalar get_emissivity_abe_matsui( const AtmosphereParameters *, Atmosphere * );
PetscScalar get_emissivity_from_flux( const Atmosphere *, const AtmosphereParameters *, PetscScalar );
PetscErrorCode set_atmosphere_volatile_content( const AtmosphereParameters *, Atmosphere * );
PetscErrorCode JSON_add_atmosphere( DM dm, const Parameters *, const Atmosphere *, const char *, cJSON *);
PetscScalar get_initial_volatile( const AtmosphereParameters *Ap, const VolatileParameters * );
PetscScalar get_dxdt( const AtmosphereParameters *Ap, const Atmosphere *, const VolatileParameters *, const Volatile * );

#endif
