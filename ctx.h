#ifndef CTX_H_
#define CTX_H_

#include <petsc.h>
#include "parameters.h"
#include "dimensionalisablefield.h"
#include "rheologicalfront.h"
#include "atmosphere.h"

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

#define NUMSOLUTIONVECS_B 43
#define NUMSOLUTIONVECS_S 30
typedef struct Solution_ {

    DimensionalisableField solutionFields_b[NUMSOLUTIONVECS_B];
    DimensionalisableField solutionFields_s[NUMSOLUTIONVECS_S];

    // TODO: eventually get rid of these
    Vec alpha, alpha_mix, cond, cp, cp_mix, dfusdr, dfusdr_temp, dSdr, dSliqdr, dSsoldr, dTdrs, dTdrs_mix, Etot, fusion, fusion_curve, fusion_curve_temp, fusion_rho, fusion_temp, fwtl, fwts, gphi, gsuper, Jcond, Jconv, Jgrav, Jmix, Jtot, kappac, kappah, liquidus, liquidus_rho, liquidus_temp, nu, phi, Ra, regime, rho, S, solidus, solidus_rho, solidus_temp, temp, visc;

    // TODO: eventually get rid of these
    Vec cp_s, cp_mix_s, dSdt_s, fusion_s, fusion_curve_s, fusion_curve_temp_s, fusion_temp_s, fwtl_s, fwts_s, gphi_s, Hradio_s, Htidal_s, Htot_s, capacitance_s, liquidus_rho_s, liquidus_s, liquidus_temp_s, phi_s, rho_s, S_s, solidus_s, solidus_rho_s, solidus_temp_s, temp_s, Hal26_s, Hk40_s, Hfe60_s, Hth232_s, Hu235_s, Hu238_s;

} Solution;

/* Some helpers for keeping track of sub-fields
 - To add a new field type, update these three things!
 - Note that these need to start with 0 to match the arrays underneath */
static const PetscInt SPIDER_NUM_FIELD_IDS = 4;
typedef enum {
  SPIDER_SOLUTION_FIELD_UNDEFINED     = 0,
  SPIDER_SOLUTION_FIELD_DSDR_B        = 1,
  SPIDER_SOLUTION_FIELD_S0            = 2,
  SPIDER_SOLUTION_FIELD_MO_VOLATILES  = 3,
  SPIDER_SOLUTION_FIELD_MO_REACTIONS  = 4,
} SpiderSolutionFieldID;
static const char * const SpiderSolutionFieldDescriptions[] = { "Undefined! Error!", "dS/dr","S at surface","Volatile partial pressure","Reaction total mass"}; /* Order must match the enum! */
static const char * const SpiderSolutionFieldUnits[]        = { "Undefined! Error!", "J kg$^{-1}$ K$^{-1}$ m$^{-1}$", "J kg$^{-1}$ K$^{-1}$", "Pa", "kg"}; /* Order must match the enum! */

/* A (temporary) struct that is used to set the eos properties at a
   given V,T or P, S.  There should be as  */
typedef struct EosEval_ {
  PetscScalar P; /* pressure */
  PetscScalar S; /* entropy */
  PetscScalar V; /* volume */
  PetscScalar T; /* temperature */
  PetscScalar Cp; /* heat capacity at constant pressure */
  PetscScalar Cv; /* heat capacity at constant volume */
  PetscScalar alpha; /* thermal expansion */
  PetscScalar rho; /* density */
  PetscScalar dTdPs; /* adiabatic temperature gradient */
  /* TODO could consider adding the following, although these are
     not usually derived from an EOS */
  /* PetscScalar cond; */
  /* PetscScalar visc */
} EosEval;

/* A Context for the Solver */
typedef struct Ctx_ {
  Mesh                   mesh;
  Solution               solution;
  DM                     da_b,da_s,da_point,da_volatiles,da_reactions; /* DMs for different subdomains */
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

  /* TODO: check with PS if this is the best way */
  /* this is a struct to enable me to pass current P, S conditions and
     update all material propoerties, according to a chosen EOS model.
     Currently need one for solid and one for melt, but could extend to
     multiple phases */
  EosEval                eos1_eval; /* melt */
  EosEval                eos2_eval; /* solid */

  /* "local" work vectors */
  Vec work_local_s,work_local_b;

  /* Pointer for PostStep reference data */
  void *postStepData;

  /* Control flags */
  PetscBool stopEarly;

  // TODO: there is no real reason that we have some data here and some in Solution and some in Mesh.
  //       it could all be flattened out (and probably would be in a re-implementation)
} Ctx;
// TODO: would be cleaner if Ctx were a pointer to _p_Ctx, as with PETSc objects.

PetscErrorCode SetupCtx(Ctx *);
PetscErrorCode DestroyCtx(Ctx *);

#endif
