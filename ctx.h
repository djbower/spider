#ifndef CTX_H_
#define CTX_H_

#include <petsc.h>
#include "parameters.h"
#include "dimensionalisablefield.h"
#include "rheologicalfront.h"
#include "atmosphere.h"
#include "eos.h"

/* common structures */

#define NUMMESHVECS_B 7
#define NUMMESHVECS_S 7
typedef struct Mesh_ {

    DimensionalisableField meshFields_b[NUMMESHVECS_B];
    DimensionalisableField meshFields_s[NUMMESHVECS_S];

    // eventually get rid of these (though think about how to better name the DimensionalisableFields..)
    Vec area_b,dPdr_b,pressure_b,radius_b,mix_b,layer_b,xi_b,dxidr_b;
    Vec pressure_s,radius_s,volume_s,dPdr_s,area_s,mass_s,xi_s;

    PetscScalar mantle_mass;

} Mesh;

#define NUMSOLUTIONVECS_B 22
#define NUMSOLUTIONVECS_S 10
typedef struct Solution_ {

    DimensionalisableField solutionFields_b[NUMSOLUTIONVECS_B];
    DimensionalisableField solutionFields_s[NUMSOLUTIONVECS_S];

    // TODO: eventually get rid of these
    Vec alpha, cond, cp, dSdxi, dTdxis, Etot, gsuper, Jcond, Jconv, Jgrav, Jmix, Jtot, kappac, kappah, nu, phi, Ra, regime, rho, S, temp, visc;

    // TODO: eventually get rid of these
    Vec cp_s, dSdt_s, Hradio_s, Htidal_s, Htot_s, capacitance_s, phi_s, rho_s, S_s, temp_s;

} Solution;

/* Some helpers for keeping track of sub-fields
 - To add a new field type, update these three things!
 - Note that these need to start with 0 to match the arrays underneath */
static const PetscInt SPIDER_NUM_FIELD_IDS = 4;
typedef enum {
  SPIDER_SOLUTION_FIELD_UNDEFINED     = 0,
  SPIDER_SOLUTION_FIELD_DSDXI_B       = 1,
  SPIDER_SOLUTION_FIELD_S0            = 2,
  SPIDER_SOLUTION_FIELD_MO_VOLATILES  = 3,
  SPIDER_SOLUTION_FIELD_MO_REACTIONS  = 4,
} SpiderSolutionFieldID;
static const char * const SpiderSolutionFieldDescriptions[] = { "Undefined! Error!", "dS/dxi","S at top staggered node","Volatile partial pressure","Reaction mass"}; /* Order must match the enum! */
static const char * const SpiderSolutionFieldUnits[]        = { "Undefined! Error!", "J kg$^{-1}$ K$^{-1}$ m$^{-1}$", "J kg$^{-1}$ K$^{-1}$", "Pa", "kg"}; /* Order must match the enum! */
static const char * const SpiderRhsFieldUnits[]             = { "Undefined! Error!", "J kg$^{-1}$ K$^{-1}$ m$^{-1}$ s$^{-1}$", "J kg$^{-1}$ K$^{-1}$ s$^{-1}$", "Pa s$^{-1}$", "kg s$^{-1}$"}; /* Order must match the enum! */

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
  Parameters             parameters;
  DimensionalisableField solDF, rhsDF; /* The solution and attached scalings */
  RheologicalFront       rheological_front_phi;
  RheologicalFront       rheological_front_dynamic;

  /* "local" work vectors */
  Vec work_local_s,work_local_b;

  /* for passing current solver time to objective_function_current_state */
  PetscReal t;

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
