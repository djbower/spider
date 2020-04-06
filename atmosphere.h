#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "parameters.h"
#include "cJSON.h"
#include "dimensionalisablefield.h"

// Note: the practice here of mixing state and parameters is suboptimal design and shouldn't be imitated.

/* these structures hold calculated quantities */
/* the parameter equivalents are in parameters.h */

typedef struct Volatile_ {
    PetscScalar x; // ppm in liquid mantle
    PetscScalar p; // partial pressure at surface (Pa)
    PetscScalar dpdt;
    PetscScalar dxdp; // dx/dp (mass fraction/Pa)
    PetscScalar mass_atmos; // mass in atmosphere (kg)
    PetscScalar mass_liquid; // mass in liquid (kg)
    PetscScalar mass_solid; // mass in solid (kg)
    PetscScalar mass_reaction; // mass exchange (gain or loss) due to reactions (kg)
    PetscScalar tau; // optical_depth at surface (non-dimensional)
    PetscScalar mixing_ratio;
    PetscScalar column_density;
    PetscScalar Knudsen; // Knudsen number
    PetscScalar jeans; // surface Jeans parameter
    PetscScalar f_thermal_escape;
    PetscScalar R_thermal_escape;
} Volatile;

typedef struct Reaction_ {
    PetscScalar dmrdt;
} Reaction;

#define NUMATMSTRUCTVECS 4
typedef struct Atmosphere_ {
    /* TODO: some of these quantities are not really strictly
       related to the atmosphere, and should perhaps live
       elsewhere */
    PetscScalar Mliq; // mass of liquid (kg)
    PetscScalar Msol; // mass of solid (kg)
    PetscScalar dMliqdt; // dMliq/dt (kg/yr)
    PetscScalar tsurf; // surface temperature
    PetscScalar dtsurfdt; // time derivative of surface temperature
    PetscScalar log10fO2; // oxygen fugacity (non-dimensional)
    PetscScalar dlog10fO2dT; // temp derivative of fO2
    PetscScalar psurf; // surface pressure
    PetscScalar dpsurfdt; // time derivative of surface pressure
    PetscScalar tau; // aggregate optical depth at surface (non-dimensional)
    PetscScalar Fatm; // net upward atmosphere flux
    PetscScalar emissivity; // variable emissivity (see also EMISSIVITY0 in AtmosphereParameters)
    Volatile    volatiles[SPIDER_MAX_VOLATILE_SPECIES]; // volatile quantities
    PetscScalar molar_mass; // mean molar mass
    Reaction    reactions[SPIDER_MAX_REACTIONS]; // reaction quantities
    PetscScalar mass_reaction[SPIDER_MAX_REACTIONS];
    DM          da_atm; // da for outputing atmosphere structure (below)
    DimensionalisableField atm_struct[NUMATMSTRUCTVECS];
    Vec atm_struct_tau;
    Vec atm_struct_temp;
    Vec atm_struct_pressure;
    Vec atm_struct_depth;
} Atmosphere;

PetscErrorCode initialise_atmosphere( Atmosphere *, const AtmosphereParameters *, const ScalingConstants );
PetscErrorCode destroy_atmosphere( Atmosphere * );

PetscScalar get_grey_body_flux( const Atmosphere *, const AtmosphereParameters * );
PetscScalar get_steam_atmosphere_zahnle_1988_flux( const Atmosphere *, const ScalingConstants );
PetscScalar get_emissivity_abe_matsui( Atmosphere *, const AtmosphereParameters *);
PetscScalar get_residual_volatile_mass( Atmosphere *, const AtmosphereParameters *, const VolatileParameters *,  const Volatile *);
PetscScalar get_emissivity_from_flux( const Atmosphere *, const AtmosphereParameters *, PetscScalar );
PetscErrorCode set_surface_temperature_from_flux( Atmosphere *, const AtmosphereParameters * );
PetscErrorCode set_reservoir_volatile_content( Atmosphere *, const AtmosphereParameters * );
PetscErrorCode set_volatile_abundances_from_partial_pressure( Atmosphere *, const AtmosphereParameters * );
PetscErrorCode JSON_add_atmosphere( DM dm, const Parameters, Atmosphere *, const char *, cJSON *);
PetscErrorCode objective_function_volatile_evolution( SNES, Vec, Vec, void * );
PetscScalar get_dpdt( Atmosphere *, const AtmosphereParameters *, PetscInt, const PetscScalar * );
PetscErrorCode set_oxygen_fugacity( Atmosphere *, const AtmosphereParameters *, const ScalingConstants );

#endif
