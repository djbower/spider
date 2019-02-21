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
    PetscScalar p; // partial pressure at surface (Pa)
    PetscScalar dpdx; // dp/dx (Pa/mass fraction)
    PetscScalar m; // mass in atmosphere (kg)
    PetscScalar tau; // optical_depth at surface (non-dimensional)
    PetscScalar mixing_ratio;
} Volatile;

#define NUMATMSTRUCTVECS 4
typedef struct Atmosphere_ {
    PetscScalar Mliq; // mass of liquid (kg)
    PetscScalar Msol; // mass of solid (kg)
    PetscScalar dMliqdt; // dMliq/dt (kg/yr)
    PetscScalar tsurf; // surface temperature
    PetscScalar tau; // aggregate optical depth at surface (dimensionless)
    PetscScalar Fatm; // net upward atmosphere flux
    PetscScalar emissivity; // variable emissivity (see also EMISSIVITY0 in AtmosphereParameters)
    Volatile    CO2; // CO2 volatile quantities
    Volatile    H2O; // H2O volatile quantities
    PetscScalar molecular_mass; // mean molecular mass
    DM          da_atm; // da for outputing atmosphere structure (below)
    DimensionalisableField atm_struct[NUMATMSTRUCTVECS];
    Vec atm_struct_tau;
    Vec atm_struct_temp;
    Vec atm_struct_pressure;
    Vec atm_struct_depth;
} Atmosphere;

PetscErrorCode initialise_atmosphere( Atmosphere *, const Constants *);
PetscErrorCode destroy_atmosphere( Atmosphere * );

PetscScalar get_grey_body_flux( const Atmosphere *, const AtmosphereParameters * );
PetscScalar get_steam_atmosphere_zahnle_1988_flux( const Atmosphere *, const Constants *C );
PetscScalar get_emissivity_abe_matsui( const AtmosphereParameters *, Atmosphere * );
#if 0
// no longer used, but kept in case it is needed in the future
PetscScalar get_emissivity_from_flux( const Atmosphere *, const AtmosphereParameters *, PetscScalar );
#endif
PetscErrorCode set_surface_temperature_from_flux( Atmosphere *, const AtmosphereParameters * );
PetscErrorCode set_atmosphere_volatile_content( const AtmosphereParameters *, Atmosphere * );
PetscErrorCode JSON_add_atmosphere( DM dm, const Parameters *, Atmosphere *, const char *, cJSON *);
// FIXME: needs replacing with PETSc non-linear solver
//PetscScalar get_initial_volatile( const AtmosphereParameters *Ap, const VolatileParameters * );
PetscScalar get_dxdt( const AtmosphereParameters *Ap, const Atmosphere *, const VolatileParameters *, const Volatile * );

#endif
