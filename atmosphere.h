#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include <petsc.h>
#include "parameters.h"
#include "cJSON.h"
#include "dimensionalisablefield.h"

/* these structures hold calculated quantities */
/* the parameter equivalents are in parameters.h */

typedef struct Volatile_ {
    PetscScalar x; // ppm in liquid mantle
    PetscScalar p; // partial pressure (Pa)
    PetscScalar dpdx; // dp/dx (Pa/mass fraction)
    PetscScalar m; // mass in atmosphere (kg)
    PetscScalar tau; // optical_depth (non-dimensional)
} Volatile;

typedef struct Atmosphere_ {
    PetscScalar Mliq; // mass of liquid (kg)
    PetscScalar Msol; // mass of solid (kg)
    PetscScalar dMliqdt; // dMliq/dt (kg/yr)
    PetscScalar tsurf; // surface temperature
    PetscScalar tau; // aggregate optical depth (dimensionless)
    PetscScalar Fatm; // net upward atmosphere flux
    PetscScalar emissivity; // variable emissivity (see also EMISSIVITY0 in AtmosphereParameters)
    Volatile    CO2; // CO2 volatile quantities
    Volatile    H2O; // H2O volatile quantities
    DM const *  da_atm_ptr;
    DimensionalisableField atm_struct[4];
    Vec atm_struct_tau;
    Vec atm_struct_temp;
    Vec atm_struct_pressure;
    Vec atm_struct_depth;
} Atmosphere;

PetscErrorCode initialise_atmosphere( DM const *, Atmosphere *, const Constants *);
PetscErrorCode destroy_atmosphere( Atmosphere * );

PetscScalar get_grey_body_flux( const Atmosphere *, const AtmosphereParameters * );
PetscScalar get_steam_atmosphere_zahnle_1988_flux( const Atmosphere *, const Constants *C );
PetscScalar get_emissivity_abe_matsui( const AtmosphereParameters *, Atmosphere * );
PetscErrorCode set_atm_struct( const AtmosphereParameters *, Atmosphere * );
PetscScalar get_emissivity_from_flux( const Atmosphere *, const AtmosphereParameters *, PetscScalar );
PetscErrorCode set_atmosphere_volatile_content( const AtmosphereParameters *, Atmosphere * );
PetscErrorCode JSON_add_atmosphere( DM dm, const Parameters *, const Atmosphere *, const char *, cJSON *);
PetscScalar get_initial_volatile( const AtmosphereParameters *Ap, const VolatileParameters * );
PetscScalar get_dxdt( const AtmosphereParameters *Ap, const Atmosphere *, const VolatileParameters *, const Volatile * );

#endif
