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


PetscErrorCode set_initial_xCO2( Atmosphere *, AtmosphereParameters const * );
PetscErrorCode set_initial_xH2O( Atmosphere *, AtmosphereParameters const * );
PetscErrorCode set_emissivity_abe_matsui( Atmosphere *, AtmosphereParameters const * );
PetscErrorCode set_dx0dt( Atmosphere *, AtmosphereParameters const *);
PetscErrorCode set_dx1dt( Atmosphere *, AtmosphereParameters const *);

#endif
