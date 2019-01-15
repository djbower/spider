#ifndef CTX_H_
#define CTX_H_

#include "petsc.h"
#include "parameters.h"
#include "dimensionalisablefield.h"

/* common structures */

#define NUMMESHVECS_B 6
#define NUMMESHVECS_S 7
typedef struct Mesh_ {

    DimensionalisableField meshFields_b[NUMMESHVECS_B];
    DimensionalisableField meshFields_s[NUMMESHVECS_S];

    // TODO: eventually get rid of these (though think about how to better name the DimensionalisableFields..)
    Vec area_b,dPdr_b,pressure_b,radius_b,mix_b,layer_b;

    // TODO: eventually get rid of these
    Vec pressure_s,radius_s,volume_s,dPdr_s,area_s,rho_s,mass_s;

    /* DJB atmosphere.  For seeing what the 'pressure' estimate of the
       mass is */
    PetscScalar mantle_mass;

} Mesh;

#define NUMSOLUTIONVECS_B 42
#define NUMSOLUTIONVECS_S 30
typedef struct Solution_ {

    DimensionalisableField solutionFields_b[NUMSOLUTIONVECS_B];
    DimensionalisableField solutionFields_s[NUMSOLUTIONVECS_S];

    // TODO: eventually get rid of these
    Vec alpha, alpha_mix, cond, cp, cp_mix, dfusdr, dfusdr_temp, dSdr, dSliqdr, dSsoldr, dTdrs, dTdrs_mix, Etot, fusion, fusion_curve, fusion_curve_temp, fusion_rho, fusion_temp, fwtl, fwts, gphi, gsuper, Jcond, Jconv, Jgrav, Jmix, Jtot, kappac, kappah, liquidus, liquidus_rho, liquidus_temp, nu, phi, regime, rho, S, solidus, solidus_rho, solidus_temp, temp, visc;

    // TODO: eventually get rid of these
    Vec cp_s, cp_mix_s, dSdt_s, fusion_s, fusion_curve_s, fusion_curve_temp_s, fusion_temp_s, fwtl_s, fwts_s, gphi_s, Hradio_s, Htidal_s, Htot_s, lhs_s, liquidus_rho_s, liquidus_s, liquidus_temp_s, phi_s, rho_s, S_s, solidus_s, solidus_rho_s, solidus_temp_s, temp_s, Hal26_s, Hk40_s, Hfe60_s, Hth232_s, Hu235_s, Hu238_s;

} Solution;

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

/* rheological front */
typedef struct RheologicalFrontMantleProperties_ {
    PetscScalar phi;
    PetscScalar depth;
    PetscScalar pressure;
    PetscScalar temperature;
} RheologicalFrontMantleProperties;

typedef struct RheologicalFront_ {
    PetscScalar phi_critical; // user-defined
    PetscInt mesh_index; // updated by code during time stepping
    PetscScalar depth; // updated by code during time stepping
    PetscScalar pressure; // updated by code during time stepping
    RheologicalFrontMantleProperties above_middle;
    RheologicalFrontMantleProperties above_mass_avg;
    RheologicalFrontMantleProperties below_middle;
    RheologicalFrontMantleProperties below_mass_avg;
} RheologicalFront;

/* Some helpers for keeping track of sub-fields
 - To add a new field type, update these three things!
*/
static const PetscInt SPIDER_NUM_FIELD_IDS=5;
typedef enum {
  SPIDER_SOLUTION_FIELD_UNDEFINED=0,
  SPIDER_SOLUTION_FIELD_DSDR_B,
  SPIDER_SOLUTION_FIELD_S0,
  SPIDER_SOLUTION_FIELD_MO_CO2,
  SPIDER_SOLUTION_FIELD_MO_H2O
} SpiderSolutionFieldID;
static const char * const SpiderSolutionFieldDescriptions[] = { "Undefined! Error!", "dS/dr (basic nodes)","S at surface","Magma Ocean CO2 content","Magma Ocean H2O content"}; /* Order must match the enum! */
static const char * const SpiderSolutionFieldUnits[] = { "Undefined! Error!", "J kg$^{-1}$ K$^{-1}$ m$^{-1}$", "J kg$^{-1}$ K$^{-1}$", "ppm", "ppm"}; // Order must match the enum! */

/* A Context for the Solver */
typedef struct Ctx_ {
  Mesh                   mesh;
  Solution               solution;
  DM                     da_b,da_s,da_point; /* DMs for different subdomains */
  DM                     dm_sol; /* A composite DM for all fields used in the timestepper */
  PetscInt               numFields; /* Number of sub-DMs in dm_sol */
  SpiderSolutionFieldID  *solutionFieldIDs; /* Which fields are in which slot */
  PetscInt               *solutionSlots;    /* The inverse map */
  Atmosphere             atmosphere;
  Mat                    qty_at_b, ddr_at_b;
  Parameters             parameters;
  DimensionalisableField solDF; /* The solution and attached scalings */
  RheologicalFront       rheological_front_phi;
  RheologicalFront       rheological_front_dynamic;

  /* "local" work vectors */
  Vec work_local_s,work_local_b;

  // TODO: there is no real reason that we have some data here and some in Solution and some in Mesh.
  //       it could all be flattened out (and probably would be in a re-implementation)
} Ctx;
// TODO: would be cleaner if Ctx were a pointer to _p_Ctx, as with PETSc objects.

PetscErrorCode SetupCtx(Ctx *);
PetscErrorCode DestroyCtx(Ctx *);

#endif
