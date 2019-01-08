#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <petsc.h>
#include "lookup.h"

/* dimensionalising constants */
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

/* atmosphere */
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

/* for storing atmosphere outputs */
typedef enum {MO_ATMOSPHERE_TYPE_GREY_BODY=1,MO_ATMOSPHERE_TYPE_ZAHNLE,MO_ATMOSPHERE_TYPE_VOLATILES,MO_ATMOSPHERE_TYPE_HEAT_FLUX,MO_ATMOSPHERE_TYPE_ENTROPY} MagmaOceanAtmosphereType;
typedef struct AtmosphereParameters_ {
    // input parameters
    MagmaOceanAtmosphereType SURFACE_BC;
    PetscBool HYBRID;
    PetscScalar surface_bc_value;
    // below are standard, also used for grey-body atmosphere
    PetscScalar emissivity0;
    PetscScalar sigma;
    PetscScalar teqm;
    PetscBool   PARAM_UTBL;
    PetscScalar param_utbl_const;
    // for volatile ODE
    PetscBool SOLVE_FOR_VOLATILES;
    PetscScalar P0; 
    VolatileParameters H2O_volatile_parameters;
    VolatileParameters CO2_volatile_parameters;
} AtmosphereParameters;

/* radiogenic heating */
typedef struct RadiogenicIsotopeParameters_ {
    PetscScalar t0;
    PetscScalar abundance; // isotopic abundance at t0
    PetscScalar concentration; // elemental concentration at t0
    PetscScalar heat_production;
    PetscScalar half_life;
} RadiogenicIsotopeParameters;

/* compositional differentiation */
typedef struct CompositionalParameters_ {
    PetscScalar X0Brg; // user-defined
    PetscScalar mass_ratio_liquidus; // computed by code
    PetscScalar muRes_muBrg; // user-defined
    /* next should really be the same as phi_critical, unless you
       have a really convincing reason otherwise! */
    PetscScalar rheological_front_phi; // user-defined
    PetscInt rheological_front_index; // updated by code during time stepping
    PetscScalar rheological_front_depth; // updated by code during time stepping
    PetscScalar rheological_front_pressure; // updated by code during time stepping
    PetscScalar mo_crystal_fraction; // updated by code during time stepping
    PetscScalar mo_bridgmanite_fraction; // updated by code during time stepping
    PetscScalar mo_mass_ratio; // updated by code during time stepping
} CompositionalParameters;

typedef enum {MO_CORE_TYPE_COOLING=1,MO_CORE_TYPE_HEAT_FLUX,MO_CORE_TYPE_ENTROPY} MagmaOceanCoreType;
typedef struct _Parameters {

    // Discretization parameters
    PetscInt    nstepsmacro,maxsteps;
    PetscReal   dtmacro;
    PetscReal   dtmacro_years; // required for output
    PetscReal   t0; /* Initial time */
    PetscInt    numpts_b,numpts_s;

    // Output Parameters
    PetscBool monitor;
    char outputDirectory[PETSC_MAX_PATH_LEN];

    // Lookup tables
    char        liquidusFilename[PETSC_MAX_PATH_LEN];
    char        solidusFilename[PETSC_MAX_PATH_LEN];
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
    Lookup      melt_prop;
    Lookup      solid_prop;

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
    PetscInt initial_condition;
    char ic_filename[PETSC_MAX_PATH_LEN];
    PetscScalar ic_melt_pressure;
    PetscScalar ic_adiabat_entropy; // entropy at top of adiabat
    PetscScalar ic_surface_entropy; // initial entropy at surface
    PetscScalar ic_core_entropy; // initial entropy at core-mantle boundary
    PetscScalar ic_dsdr;
    PetscScalar radius;
    PetscScalar coresize;
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

    RadiogenicIsotopeParameters al26_parameters;
    RadiogenicIsotopeParameters k40_parameters;
    RadiogenicIsotopeParameters fe60_parameters;
    RadiogenicIsotopeParameters th232_parameters;
    RadiogenicIsotopeParameters u235_parameters;
    RadiogenicIsotopeParameters u238_parameters;

    // Composition
    CompositionalParameters compositional_parameters;

} Parameters;

PetscErrorCode InitializeParametersAndSetFromOptions(Parameters *parameters);
PetscErrorCode PrintParameters(Parameters const *parameters);
PetscErrorCode SetLookups( Parameters * );

#endif
