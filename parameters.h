#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <petsc.h>

/*
 ******************************************************************************
 * Dimensional constants
 ******************************************************************************
 */

typedef struct _Constants {
    // primary
    PetscScalar RADIUS;
    PetscScalar TEMP;
    PetscScalar ENTROPY;
    PetscScalar DENSITY;
    PetscScalar VOLATILE;
    // derived from primary
    PetscScalar AREA;
    PetscScalar VOLUME;
    PetscScalar MASS;
    PetscScalar TIME;
    PetscScalar TIMEYRS;
    PetscScalar SENERGY;
    PetscScalar ENERGY;
    PetscScalar PRESSURE;
    PetscScalar POWER;
    PetscScalar FLUX;
    PetscScalar DPDR;
    PetscScalar GRAVITY;
    PetscScalar KAPPA;
    PetscScalar DTDP;
    PetscScalar DSDR;
    PetscScalar DTDR;
    PetscScalar GSUPER;
    PetscScalar VISC;
    PetscScalar LOG10VISC;
    PetscScalar COND;
    PetscScalar SIGMA;
    PetscScalar LHS;
    PetscScalar RHS;
    PetscScalar HEATGEN;
} Constants;

/*
 ******************************************************************************
 * Equation of state parameters
 ******************************************************************************
 */

/* 1-D lookup */
typedef struct _Interp1d {
    PetscInt    NX; 
    PetscScalar *xa;
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar *ya;
    PetscScalar ymin;
    PetscScalar ymax;
} Interp1d;

/* 2-D lookup */
typedef struct _Interp2d {
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
} Interp2d;

/* lookup */
typedef struct _Lookup {
    /* lookup data filenames */
    char        rho_filename[PETSC_MAX_PATH_LEN];
    char        dTdPs_filename[PETSC_MAX_PATH_LEN];
    char        cp_filename[PETSC_MAX_PATH_LEN];
    char        temp_filename[PETSC_MAX_PATH_LEN];
    char        alpha_filename[PETSC_MAX_PATH_LEN];
    /* TODO: eventually melting curves may move elsewhere */
    char        liquidus_filename[PETSC_MAX_PATH_LEN];
    char        solidus_filename[PETSC_MAX_PATH_LEN];
    /* lookup objects to evaluate */
    Interp2d rho; /* density, kg/m^3 */
    Interp2d dTdPs; /* adiabatic temperature gradient, K/Pa */
    Interp2d cp; /* heat capacity, J/kg/K */
    Interp2d temp; /* temperature, K */
    Interp2d alpha; /* thermal expansion, 1/K */
    PetscScalar cond; /* thermal conductivity, W/m/K */
    PetscScalar log10visc; /* log base 10 of viscosity */
    /* TODO: might be better to store melting curve data
       in a different struct */
    Interp1d liquidus; /* liquidus, J/kg/K */
    Interp1d solidus; /* solidus, J/kg/K */
} Lookup;

/* EOS */
typedef struct _EosParameters {
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
    PetscScalar Avogadro;
    PetscScalar PV_UNIT;
    PetscScalar KBOLTZ;
    PetscScalar bscale;
    Lookup lookup;
} EosParameters;

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
typedef struct VolatileParameters_ {
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
} VolatileParameters;

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
typedef struct AtmosphereParameters_ {
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
    PetscScalar kB;
    PetscScalar bigG;
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
    PetscScalar         Rgas; // gas constant
    PetscScalar const * gravity_ptr;
    PetscScalar const * radius_ptr;
    PetscScalar const * VOLATILE_ptr;
    PetscScalar const * mantle_mass_ptr;
} AtmosphereParameters;

/*
 ******************************************************************************
 * Radiogenic heating parameters
 ******************************************************************************
 */

typedef struct RadiogenicIsotopeParameters_ {
    PetscScalar t0;
    PetscScalar abundance; // isotopic abundance at t0
    PetscScalar concentration; // elemental concentration at t0
    PetscScalar heat_production;
    PetscScalar half_life;
} RadiogenicIsotopeParameters;

#if 0
/* compositional differentiation */
typedef struct CompositionParameters_ {
    PetscScalar Brg_initial_fraction; // user-defined
    PetscScalar Res_Brg_mass_ratio; // user-defined
    /* next scales the density at the liquidus */
    PetscScalar BSE_Brg_mass_ratio_at_liquidus; // computed by code, but fixed with time
} CompositionParameters;
#endif

/*
 ******************************************************************************
 * Main parameters structure
 ******************************************************************************
*/

typedef enum {MO_CORE_TYPE_COOLING=1,MO_CORE_TYPE_HEAT_FLUX,MO_CORE_TYPE_ENTROPY} MagmaOceanCoreType;
typedef struct _Parameters {

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

    PetscInt    SOLID_EOS;
    PetscInt    MELT_EOS;

// FIXME: MOVED INTO EOS PARAMETERS, LOOKUP
    // Lookup tables
    char        liquidusFilename[PETSC_MAX_PATH_LEN];
    char        solidusFilename[PETSC_MAX_PATH_LEN];
#if 0
    char        alphaSolFilename[PETSC_MAX_PATH_LEN];
    char        alphaMelFilename[PETSC_MAX_PATH_LEN];
    char        cpSolFilename[PETSC_MAX_PATH_LEN];
    char        cpMelFilename[PETSC_MAX_PATH_LEN];
    char        dtdpsSolFilename[PETSC_MAX_PATH_LEN];
    char        dtdpsMelFilename[PETSC_MAX_PATH_LEN];
    char        rhoSolFilename[PETSC_MAX_PATH_LEN];
    char        rhoMelFilename[PETSC_MAX_PATH_LEN];
    char        tempSolFilename[PETSC_MAX_PATH_LEN];
    char        tempMelFilename[PETSC_MAX_PATH_LEN];
#endif

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
    PetscBool COMPOSITION; // Brg and Res compositional model
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
    PetscScalar log10visc_sol;
    PetscScalar activation_energy_sol;
    PetscScalar activation_volume_sol;
    PetscScalar log10visc_min;
    PetscScalar log10visc_max;
    PetscScalar Mg_Si0;
    PetscScalar Mg_Si1;
    PetscScalar cond_sol;
    PetscScalar log10visc_mel;
    PetscScalar cond_mel;
    PetscScalar eddy_diffusivity_thermal;
    PetscScalar eddy_diffusivity_chemical;
    PetscInt    VISCOUS_LID;
    PetscScalar lid_log10visc;
    PetscScalar lid_thickness;

    // Additional Atmosphere Parameters
    AtmosphereParameters atmosphere_parameters;

    // Scaling factors / dimensionalising constants
    Constants constants;

    // Radiogenic heating parameters
    RadiogenicIsotopeParameters al26_parameters;
    RadiogenicIsotopeParameters k40_parameters;
    RadiogenicIsotopeParameters fe60_parameters;
    RadiogenicIsotopeParameters th232_parameters;
    RadiogenicIsotopeParameters u235_parameters;
    RadiogenicIsotopeParameters u238_parameters;

    // equation of state for melt and solid
    EosParameters eos1_parameters; /* for melt */
    EosParameters eos2_parameters; /* for solid */

    // FIXME
    // Composition
    //CompositionParameters composition_parameters;

} Parameters;

/* Parameters Methods */
PetscErrorCode InitializeParametersAndSetFromOptions( Parameters * );
PetscErrorCode PrintParameters( Parameters const * );
PetscErrorCode ParametersDestroy( Parameters * );

/* ReactionParameters Methods */
PetscErrorCode ReactionParametersCreateAmmonia1( ReactionParameters *, const AtmosphereParameters * );
PetscErrorCode ReactionParametersCreateCarbonDioxide1( ReactionParameters *, const AtmosphereParameters * );
PetscErrorCode ReactionParametersCreateMethane1( ReactionParameters *, const AtmosphereParameters * );
PetscErrorCode ReactionParametersCreateWater1( ReactionParameters *, const AtmosphereParameters * );
PetscErrorCode ReactionParametersCreateSimpleWater1( ReactionParameters *, const AtmosphereParameters * );
PetscErrorCode ReactionParametersCreateSimpleWater2( ReactionParameters *, const AtmosphereParameters * );
PetscErrorCode ReactionParametersCreateSimpleWater3( ReactionParameters *, const AtmosphereParameters * );
PetscErrorCode ReactionParametersCreateMiddleWater2( ReactionParameters *, const AtmosphereParameters * );
PetscErrorCode ReactionParametersDestroy( ReactionParameters *);

#endif
