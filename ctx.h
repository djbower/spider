#ifndef CTX_H_
#define CTX_H_

#include "petsc.h"
#include "parameters.h"
#include "dimensionalisablefield.h"

/* common structures */

#define NUMMESHVECS_B 5
#define NUMMESHVECS_S 7
typedef struct Mesh_ {

    DimensionalisableField meshFields_b[NUMMESHVECS_B];
    DimensionalisableField meshFields_s[NUMMESHVECS_S];

    // PDS TODO: eventually get rid of these (though think about how to better name the DimensionalisableFields..)
    Vec area_b,dPdr_b,pressure_b,radius_b, mix_b;

    // PDS TODO: eventually get rid of these
    Vec pressure_s,radius_s,volume_s,dPdr_s,area_s,rho_s,mass_s;

    /* DJB atmosphere.  For seeing what the 'pressure' estimate of the
       mass is */
    PetscScalar mantle_mass;

} Mesh;

#define NUMSOLUTIONVECS_B 41
#define NUMSOLUTIONVECS_S 24
typedef struct Solution_ {

    DimensionalisableField solutionFields_b[NUMSOLUTIONVECS_B];
    DimensionalisableField solutionFields_s[NUMSOLUTIONVECS_S];

    // PDS TODO: eventually get rid of these
    Vec alpha, alpha_mix, cond, cp, cp_mix, dfusdr, dfusdr_temp, dSdr, dSliqdr, dSsoldr, dTdrs, dTdrs_mix, Etot, fusion, fusion_curve, fusion_curve_temp, fusion_rho, fusion_temp, fwtl, fwts, gphi, gsuper, Jcond, Jconv, Jgrav, Jmix, Jtot, kappah, liquidus, liquidus_rho, liquidus_temp, nu, phi, regime, rho, S, solidus, solidus_rho, solidus_temp, temp, visc;

    // PDS TODO: eventually get rid of these
    Vec cp_s, cp_mix_s, dSdt_s, fusion_s, fusion_curve_s, fusion_curve_temp_s, fusion_temp_s, fwtl_s, fwts_s, gphi_s, Hradio_s, Htidal_s, Htot_s, lhs_s, liquidus_rho_s, liquidus_s, liquidus_temp_s, phi_s, rho_s, S_s, solidus_s, solidus_rho_s, solidus_temp_s, temp_s;

} Solution;


/* A Context for the Solver */
typedef struct Ctx_ {
  Mesh            mesh;
  Solution        solution;
  DM              da_b,da_s,da_surface,da_mo_co2,da_mo_h2o; /* DMs for different subdomains */
  DM              dm_sol; /* A composite DM for all fields used in the timestepper */
  PetscInt        numFields; /* Number of sub-DMs in dm_sol */
  const char      **solutionFieldDescription;
  Atmosphere      atmosphere;
  Mat             qty_at_b, ddr_at_b;
  Parameters      parameters;

  /* "local" work vectors */
  Vec work_local_s,work_local_b;

  // PDS TODO:
  // Note: there is no real reason that we have some data here and some in Solution and some in Mesh.
  //       it could all be flattened out (and probably would be in a re-implementation)
} Ctx;
// TODO: would be cleaner if Ctx where a pointer to _p_Ctx, as with PETSc objects.


PetscErrorCode SetupCtx(Ctx *);
PetscErrorCode DestroyCtx(Ctx *);

#endif
