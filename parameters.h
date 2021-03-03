#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <petsc.h>
#include "constants.h"
#include "eos.h"

/*
 ******************************************************************************
 * Volatile parameters
 ******************************************************************************
 */

/* Volatiles: each has an options_prefix, set with a particular option
   and used to look for additional options to populate the other data.
   SPIDER itself is agnostic to which set of volatiles you use or how many.
   There is a hard-coded number of "slots". If it ever became important,
   this could be changed to dynamically allocate and deallocate, everywhere
   this constant is used: */
#define SPIDER_MAX_VOLATILE_SPECIES 16
typedef struct {
    char        prefix[128];  /* Maximum prefix length */
    /* next three (initial_) are not mutually exclusive */
    PetscScalar initial_total_abundance;
    PetscScalar initial_atmos_pressure;
    PetscScalar initial_ocean_moles;
    PetscScalar kdist;
    PetscScalar kabs; // note this is without pressure-dependence
    PetscInt    SOLUBILITY;
    PetscScalar henry;
    PetscScalar henry_pow;
    PetscScalar henry2;
    PetscScalar henry_pow2;
    PetscScalar jeans_value; // for thermal escape
    PetscScalar R_thermal_escape_value; // for thermal (Jean's) escape
    PetscScalar constant_escape_value; // for constant escape
    PetscScalar molar_mass;
    PetscScalar cross_section;
    PetscReal   poststep_change; // allowable fractional change (only for -activate_poststep)
} data_VolatileParameters;
typedef data_VolatileParameters* VolatileParameters;

/*
 ******************************************************************************
 * Reaction parameters
 ******************************************************************************
 */

/*
ReactionParameters: an object which describes a particular type of chemical equilibrium.

Specifically, it represents the situation where two or more volatiles can be 
converted (reversibly, unconditionally) to one another by means of chemical 
reactions, and where these volatiles are in (instantaneous, quasi-static) 
equilibrium, as described by some (time-and-state-dependent) function.

Thus, the data are:
- A set of N volatiles, identified by indices corresponding to volatile species.
- A set of N coefficients stoichiometry_i
  The sign of these defines the direction of the reaction; positive implies a product, negative a reactant.

Developers' Note: this class is done "properly", in that Reaction is a pointer
to a struct, not itself a struct.  Other classes in SPIDER aren't (yet) all
done this way.
*/
typedef struct {
  const char *type;
  PetscInt   n_volatiles;
  PetscInt   *volatiles; /* indices for volatiles. Populated from prefix strings during parameter processing */
  PetscReal  *stoichiometry;  /* stoichiometry */
  PetscReal  *Keq_coeffs;  /* Equilibrium constant parameters */
  PetscReal  fO2_stoichiometry; /* fO2 stoichiometry (zero turns it off, non-zero value is used) */
} data_ReactionParameters;
typedef data_ReactionParameters* ReactionParameters;

/* We hard-code the maximum number of reactions. If needbe, the array of
   ReactionParameters could be dynamically allocated and freed.*/
#define SPIDER_MAX_REACTIONS 8

/*
 ******************************************************************************
 * Atmosphere parameters
 ******************************************************************************
 */

typedef enum {MO_ATMOSPHERE_TYPE_GREY_BODY=1,MO_ATMOSPHERE_TYPE_ZAHNLE,MO_ATMOSPHERE_TYPE_VOLATILES,MO_ATMOSPHERE_TYPE_HEAT_FLUX,MO_ATMOSPHERE_TYPE_ENTROPY} MagmaOceanAtmosphereType;
typedef enum {OXYGEN_FUGACITY_NONE=0,OXYGEN_FUGACITY_CI,OXYGEN_FUGACITY_CV,OXYGEN_FUGACITY_H,OXYGEN_FUGACITY_EH,OXYGEN_FUGACITY_EUCRITE,OXYGEN_FUGACITY_FISCHER_IW_PLUS_HALF,OXYGEN_FUGACITY_ONEILL_IW_PLUS_HALF} OxygenFugacityType;
typedef struct {
    // input parameters
    PetscInt IC_ATMOSPHERE;
    char ic_atmosphere_filename[PETSC_MAX_PATH_LEN];
    MagmaOceanAtmosphereType SURFACE_BC;
    OxygenFugacityType OXYGEN_FUGACITY; // for chemical reactions
    PetscBool SURFACE_BC_ACC;
    /* below to replace with shallow ocean layer */
    //PetscBool VISCOUS_MANTLE_COOLING_RATE;
    PetscBool THERMAL_ESCAPE;
    PetscBool CONSTANT_ESCAPE;
    PetscScalar surface_bc_value;
    // below are standard, also used for grey-body atmosphere
    PetscScalar emissivity0;
    PetscScalar teqm;
    PetscReal   tsurf_poststep_change; // maxmimum absolute change in surface temperature before exiting loop
    PetscBool   PARAM_UTBL;
    PetscScalar param_utbl_const;
    // for volatile ODE
    PetscScalar         P0;
    PetscInt            n_volatiles;
    VolatileParameters  volatile_parameters[SPIDER_MAX_VOLATILE_SPECIES];
    PetscInt            n_reactions;
    ReactionParameters  reaction_parameters[SPIDER_MAX_REACTIONS];
    PetscScalar const * gravity_ptr;
    PetscScalar const * radius_ptr;
    PetscScalar const * VOLATILE_ptr;
    PetscScalar const * mantle_mass_ptr;
} data_AtmosphereParameters;
typedef data_AtmosphereParameters* AtmosphereParameters;

/*
 ******************************************************************************
 * Radiogenic heating parameters
 ******************************************************************************
 */

typedef struct {
    char        prefix[128];  /* Maximum prefix length */
    PetscScalar t0; // time for which concentration is known
    /* abundance and concentration are simply multiplied in the code.  This is
       because it is sometimes more convenient to give the concentration in
       ppmw of Uranium and then express the relative abundance of U235 and
       U238 */
    PetscScalar abundance; // radionuclide abundance at t0
    PetscScalar concentration; // elemental concentration at t0
    PetscScalar heat_production; // W/kg
    PetscScalar half_life; // years
} data_RadionuclideParameters;
typedef data_RadionuclideParameters* RadionuclideParameters;

/*
 ******************************************************************************
 * Main parameters structure
 ******************************************************************************
*/

#define SPIDER_MAX_RADIONUCLIDES 8
#define SPIDER_MAX_PHASES 2
typedef enum {MO_CORE_TYPE_COOLING=1,MO_CORE_TYPE_HEAT_FLUX,MO_CORE_TYPE_ENTROPY} MagmaOceanCoreType;
typedef struct {

    // Discretization parameters
    PetscInt    nstepsmacro,maxsteps;
    PetscReal   dtmacro;
    PetscReal   t0; /* initial time */
    PetscInt    numpts_b,numpts_s;

    // Rollback / Post-step logic parameters
    PetscBool   rollBackActive,postStepActive;

    // Output Parameters
    PetscBool monitor;
    char outputDirectory[PETSC_MAX_PATH_LEN];

    //  "Standard" parameters
    MagmaOceanCoreType CORE_BC;
    PetscScalar core_bc_value;

    PetscBool MASS_COORDINATES;
    PetscBool CONDUCTION;
    PetscBool CONVECTION;
    PetscBool MIXING;
    /* now an int to explore different separation schemes */
    PetscInt  SEPARATION;
    PetscBool HTIDAL;
    PetscInt mixing_length;
    PetscScalar mixing_length_a;
    PetscScalar mixing_length_b;
    PetscScalar layer_interface_radius;
    PetscInt IC_INTERIOR;
    char ic_interior_filename[PETSC_MAX_PATH_LEN];
    PetscScalar ic_melt_pressure;
    PetscScalar ic_adiabat_entropy; // entropy at top of adiabat
    PetscScalar ic_surface_entropy; // initial entropy at surface
    PetscScalar ic_core_entropy; // initial entropy at core-mantle boundary
    PetscScalar ic_dsdr;
    PetscScalar radius;
    PetscScalar coresize;
    PetscScalar coremass;
    PetscScalar rhos;
    PetscScalar beta;
    PetscScalar grain;
    PetscScalar gravity;
    PetscScalar phi_critical;
    PetscScalar phi_width;
    PetscScalar rho_core;
    PetscScalar cp_core;
    PetscScalar tfac_core_avg;
    PetscScalar activation_energy_sol;
    PetscScalar activation_volume_sol;
    PetscScalar log10visc_min;
    PetscScalar log10visc_max;
    PetscScalar Mg_Si0;
    PetscScalar Mg_Si1;
    PetscScalar eddy_diffusivity_thermal;
    PetscScalar eddy_diffusivity_chemical;
    PetscInt    VISCOUS_LID;
    PetscScalar lid_log10visc;
    PetscScalar lid_thickness;

    /* Scaling constants */
    ScalingConstants scaling_constants;

    /* Fundamental constants */
    FundamentalConstants fundamental_constants;

    /* Atmosphere parameters */
    AtmosphereParameters atmosphere_parameters;

    /* Radionuclides */
    PetscInt    n_radionuclides;
    RadionuclideParameters radionuclide_parameters[SPIDER_MAX_RADIONUCLIDES];

    /* Equation of state for mass coordinates */
    EOS      eos_mesh;

    /* Equation of state */
    EOS      eos;
    EOS      eos_phases[SPIDER_MAX_PHASES];
    PetscInt n_phases;

} data_Parameters;
typedef data_Parameters* Parameters;

/* Parameters Methods */
PetscErrorCode ParametersCreate( Parameters * );
PetscErrorCode ParametersDestroy( Parameters * );
PetscErrorCode PrintParameters( Parameters const );

/* ReactionParameters Methods */
PetscErrorCode ReactionParametersCreateAmmonia1( ReactionParameters *, const AtmosphereParameters, const ScalingConstants );
PetscErrorCode ReactionParametersCreateCarbonDioxide1( ReactionParameters *, const AtmosphereParameters, const ScalingConstants );
PetscErrorCode ReactionParametersCreateCarbonDioxideJANAF( ReactionParameters *, const AtmosphereParameters, const ScalingConstants );
PetscErrorCode ReactionParametersCreateMethane1( ReactionParameters *, const AtmosphereParameters, const ScalingConstants );
PetscErrorCode ReactionParametersCreateWaterSchaefer( ReactionParameters *, const AtmosphereParameters, const ScalingConstants );
PetscErrorCode ReactionParametersCreateWaterJANAF( ReactionParameters *, const AtmosphereParameters, const ScalingConstants );
PetscErrorCode ReactionParametersDestroy( ReactionParameters *);

#endif
