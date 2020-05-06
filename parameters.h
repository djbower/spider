#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <petsc.h>

/*
 ******************************************************************************
 * Dimensional constants
 ******************************************************************************
 */

/* constants to scale the physical problem, largely chosen based on numerical
   considerations */
typedef struct {
    /* primary */
    PetscScalar RADIUS;
    PetscScalar TEMP;
    PetscScalar ENTROPY; /* note: specific entropy */
    PetscScalar DENSITY;
    PetscScalar VOLATILE;
    /* derived from primary */
    PetscScalar AREA;
    PetscScalar VOLUME;
    PetscScalar MASS;
    PetscScalar TIME;
    PetscScalar TIMEYRS;
    PetscScalar SENERGY; /* specific energy */
    PetscScalar ENERGY;
    PetscScalar PRESSURE;
    PetscScalar POWER;
    PetscScalar FLUX;
    PetscScalar DPDR;
    PetscScalar GRAVITY;
    PetscScalar KAPPA;
    PetscScalar DSDP;
    PetscScalar DSDR;
    PetscScalar DTDP;
    PetscScalar DTDR;
    PetscScalar GSUPER;
    PetscScalar VISC;
    PetscScalar LOG10VISC;
    PetscScalar COND;
    PetscScalar SIGMA;
    PetscScalar RHS;
    PetscScalar HEATGEN;
} data_ScalingConstants;
typedef data_ScalingConstants* ScalingConstants;

/* fundamental constants */
typedef struct {
    PetscScalar AVOGADRO;
    PetscScalar BOLTZMANN;
    PetscScalar GAS;
    PetscScalar GRAVITATIONAL;
    PetscScalar STEFAN_BOLTZMANN;
} data_FundamentalConstants;
typedef data_FundamentalConstants* FundamentalConstants;

/*
 ******************************************************************************
 * Equation of state parameters
 ******************************************************************************
 */

/* 1-D lookup */
typedef struct {
    PetscInt    NX; 
    PetscScalar *xa;
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar *ya;
    PetscScalar ymin;
    PetscScalar ymax;
} data_Interp1d;
typedef data_Interp1d* Interp1d;

/* 2-D lookup */
typedef struct {
    PetscInt    NX; 
    PetscScalar *xa;
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar dx; 
    PetscInt    NY; 
    PetscScalar *ya;
    PetscScalar ymin;
    PetscScalar ymax;
    PetscScalar dy; 
    PetscScalar **za;
} data_Interp2d;
typedef data_Interp2d* Interp2d;

/* lookup */
typedef struct {
    /* lookup data filenames */
    char        rho_filename[PETSC_MAX_PATH_LEN];
    char        dTdPs_filename[PETSC_MAX_PATH_LEN];
    char        cp_filename[PETSC_MAX_PATH_LEN];
    char        temp_filename[PETSC_MAX_PATH_LEN];
    char        alpha_filename[PETSC_MAX_PATH_LEN];
    /* lookup objects to evaluate */
    /* memory is not necessarily allocated to these lookup
       objects if they are not required */
    Interp2d rho; /* density, kg/m^3 */
    Interp2d dTdPs; /* adiabatic temperature gradient, K/Pa */
    Interp2d cp; /* heat capacity, J/kg/K */
    Interp2d temp; /* temperature, K */
    Interp2d alpha; /* thermal expansion, 1/K */
} data_Lookup;
typedef data_Lookup* Lookup;

/* RTpress (Wolf and Bower, 2018) */
typedef struct {
    PetscScalar V0;
    PetscScalar T0;
    PetscScalar S0;
    PetscScalar K0;
    PetscScalar KP0;
    PetscScalar E0;
    PetscScalar gamma0;
    PetscScalar gammaP0;
    PetscScalar m;
    PetscScalar b0;
    PetscScalar b1;
    PetscScalar b2;
    PetscScalar b3;
    PetscScalar b4;
    PetscScalar mavg;
    PetscScalar PV_UNIT;
    PetscScalar KBOLTZ;
    PetscScalar bscale;
    PetscScalar const * AVOGADRO_ptr;
} data_RTpressParameters;
typedef data_RTpressParameters* RTpressParameters;

/* EOS */
/* TODO: an optional consolidation, is to combine the EosComposite structure
   with the EosParameters structure.  This can be done if the EosParameters
   structure contains an array of pointers to EosParameters (which I think
   is legal with C) */
#define SPIDER_MAX_PHASES 2
typedef struct {
    char prefix[128];  /* Maximum prefix length */
    /* Eos TYPE
       1. Lookup
       2. Rtpress
    */
    PetscInt TYPE;
    Lookup lookup;
    RTpressParameters rtpress_parameters;
    /* include thermal and transport properties */
    PetscScalar cond; /* thermal conductivity, W/m/K */
    PetscScalar log10visc; /* log base 10 of viscosity */
    /* phase boundary which is evaluated using this EOS */
    PetscBool PHASE_BOUNDARY; /* is a phase boundary for this EOS defined? */
    char phase_boundary_filename[PETSC_MAX_PATH_LEN];
    Interp1d phase_boundary; /* pressure-entropy space, J/kg/K */
} data_EosParameters;
typedef data_EosParameters* EosParameters;

#define SPIDER_MAX_COMPOSITE_PHASES 1 
typedef struct {
    const char *prefix;  /* Maximum prefix length */
    EosParameters eos_parameters[SPIDER_MAX_PHASES];
    PetscInt n_eos;
} data_EosComposite;
typedef data_EosComposite* EosComposite;

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
    PetscScalar kdist;
    PetscScalar kabs; // note this is without pressure-dependence
    PetscInt    SOLUBILITY;
    PetscScalar henry;
    PetscScalar henry_pow;
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
typedef enum {OXYGEN_FUGACITY_NONE=0,OXYGEN_FUGACITY_CI,OXYGEN_FUGACITY_CV,OXYGEN_FUGACITY_H,OXYGEN_FUGACITY_EH,OXYGEN_FUGACITY_EUCRITE,OXYGEN_FUGACITY_IW,OXYGEN_FUGACITY_IW_MINUS_ONE,OXYGEN_FUGACITY_IW_MINUS_TWO,OXYGEN_FUGACITY_IW_MINUS_THREE} OxygenFugacityType;
typedef struct {
    // input parameters
    PetscInt IC_ATMOSPHERE;
    char ic_atmosphere_filename[PETSC_MAX_PATH_LEN];
    MagmaOceanAtmosphereType SURFACE_BC;
    OxygenFugacityType OXYGEN_FUGACITY; // for chemical reactions
    PetscBool VISCOUS_MANTLE_COOLING_RATE;
    PetscBool THERMAL_ESCAPE;
    PetscBool CONSTANT_ESCAPE;
    PetscScalar surface_bc_value;
    // below are standard, also used for grey-body atmosphere
    PetscScalar emissivity0;
    PetscScalar sigma;
    PetscScalar Avogadro;
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

    PetscBool CONDUCTION;
    PetscBool CONVECTION;
    PetscBool MIXING;
    PetscBool SEPARATION;
    PetscBool HRADIO;
    PetscBool HTIDAL;
    PetscBool SOLID_CONVECTION_ONLY; // solid convection only
    PetscBool LIQUID_CONVECTION_ONLY; // liquid convection only
    PetscInt mixing_length;
    PetscScalar mixing_length_layer_radius;
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
    PetscScalar phi_skew;
    PetscScalar rho_core;
    PetscScalar cp_core;
    PetscScalar tfac_core_avg;
    PetscScalar matprop_smooth_width;
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

    /* Phases (EOS) */
    PetscInt    n_phases;
    EosParameters eos_parameters[SPIDER_MAX_PHASES];
    PetscInt    n_composite_phases;
    EosComposite eos_composites[SPIDER_MAX_COMPOSITE_PHASES];

} data_Parameters;
typedef data_Parameters* Parameters;

/* Parameters Methods */
PetscErrorCode ParametersCreate( Parameters * );
PetscErrorCode ParametersDestroy( Parameters * );
PetscErrorCode PrintParameters( Parameters const );

/* ReactionParameters Methods */
PetscErrorCode ReactionParametersCreateAmmonia1( ReactionParameters *, const AtmosphereParameters );
PetscErrorCode ReactionParametersCreateCarbonDioxide1( ReactionParameters *, const AtmosphereParameters );
PetscErrorCode ReactionParametersCreateMethane1( ReactionParameters *, const AtmosphereParameters );
PetscErrorCode ReactionParametersCreateWater1( ReactionParameters *, const AtmosphereParameters );
PetscErrorCode ReactionParametersCreateSimpleWater1( ReactionParameters *, const AtmosphereParameters );
PetscErrorCode ReactionParametersCreateSimpleWater2( ReactionParameters *, const AtmosphereParameters );
PetscErrorCode ReactionParametersCreateSimpleWater3( ReactionParameters *, const AtmosphereParameters );
PetscErrorCode ReactionParametersCreateMiddleWater2( ReactionParameters *, const AtmosphereParameters );
PetscErrorCode ReactionParametersDestroy( ReactionParameters *);

#endif
