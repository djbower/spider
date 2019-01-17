#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include <petsc.h>
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

PetscScalar get_grey_body_flux( const Atmosphere *, const AtmosphereParameters * );
PetscScalar get_steam_atmosphere_zahnle_1988_flux( const Atmosphere *, const Constants *C );
PetscScalar get_emissivity_abe_matsui( const Parameters *, Atmosphere * );
PetscScalar get_emissivity_from_flux( const Atmosphere *, const AtmosphereParameters *, PetscScalar );
PetscScalar solve_newton_method( PetscScalar, PetscScalar, PetscScalar );
PetscErrorCode set_atmosphere_volatile_content( const Parameters *, Atmosphere *, PetscScalar, PetscScalar );

#endif
