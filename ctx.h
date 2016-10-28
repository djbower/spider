#ifndef CTX_H_
#define CTX_H_

#include "petsc.h"

// Heinous. Petsc likes include headers which define I as a complex number,
// even if you aren't aware that PETSc even supports that. 
#undef I

#include "global_defs.h"

/* common structures */

/* this structure has an x and y array size equal to NLS */
typedef struct _Interp1d {
    PetscScalar xa[NLS];
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar ya[NLS];
    PetscScalar ymin;
    PetscScalar ymax;
    PetscScalar dx;
} Interp1d;

typedef struct _Interp2d {
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar ymin;
    PetscScalar ymax;
    PetscScalar xa[NX];
    PetscScalar ya[NY];
    PetscScalar za[NX*NY];
    PetscScalar dx;
    PetscScalar dy;
} Interp2d;

/* lookup for a single phase */
typedef struct _Lookup {
    Interp2d rho; /* density, kg / m^3 */
    Interp2d dTdPs; /* adiabatic temperature gradient, K / Pa */
    Interp2d cp; /* heat capacity, J / (kg K) */
    Interp2d temp; /* temperature, K */
    Interp2d alpha; /* thermal expansion, 1/K */
    PetscScalar cond; /* thermal conductivity, W / (m K) */
    PetscScalar log10visc; /* log base 10 of viscosity */
    /*for single component system these are the same for both solid
      and melt phases, but formulating it this way might provide a
      a way forward for a multicomponent system */
    Interp1d liquidus; /* liquidus, J / (kg K) */
    Interp1d solidus; /* solidus, J / (kg K) */
} Lookup;

#define NUMMESHVECS_B 5
#define NUMMESHVECS_S 5
typedef struct _Mesh {

    Vec meshVecs_b[NUMMESHVECS_B];
    Vec area_b,dPdr_b,pressure_b,radius_b, mix_b;  

    Vec meshVecs_s[NUMMESHVECS_S];
    Vec pressure_s,radius_s,volume_s,dPdr_s,area_s;

} Mesh;

#define NUMSOLUTIONVECS_B 40
#define NUMSOLUTIONVECS_S 14
typedef struct _Solution {
    /* TI means time-independent */

    Vec solutionVecs_b[NUMSOLUTIONVECS_B];
    Vec alpha,alpha_mix,cond,cp,cp_mix,dfusdr,dfusdr_temp,dphidr,dSdr,dTdPs,dTdPs_mix,dTdrs,dTdrs_mix,Etot,fusion,fusion_curve,fusion_curve_temp,fusion_rho,fusion_temp,gsuper,Jcond,Jconv,Jgrav,Jheat,Jmass,Jmix,Jtot,kappah,liquidus,liquidus_rho,liquidus_temp,nu,phi,rho,S,solidus,solidus_rho,solidus_temp,temp,visc;

    Vec solutionVecs_s[NUMSOLUTIONVECS_S];
    Vec fusion_s, fusion_curve_s, fusion_curve_temp_s, fusion_temp_s, lhs_s, liquidus_rho_s, liquidus_s, liquidus_temp_s, phi_s, rho_s, solidus_s, solidus_rho_s, solidus_temp_s, temp_s;

} Solution;

/* A Context for the Solver */
typedef struct _Ctx {
  Lookup   melt_prop;
  Lookup   solid_prop;
  Mesh     mesh;
  Solution solution;
  DM       da_b,da_s;
  PetscScalar BC_BOT_FAC; // core-cooling boundary condition
  PetscScalar S_init; // initial entropy
  Mat      qty_at_b, ddr_at_b;

  /* "local" work vectors */
  Vec work_local_s,work_local_b;

  // Note: there is no real reason that we have some data here and some in Solution and some in Mesh.
  //       it could all be flattened out (and probably would be in a re-implementation)
} Ctx;

PetscErrorCode setup_ctx(Ctx *);
PetscErrorCode destroy_ctx(Ctx *);

PetscErrorCode set_d_dr( Ctx * );
PetscErrorCode set_capacitance( Ctx *, Vec );
PetscErrorCode set_lookups( Ctx * );
PetscErrorCode set_matprop_and_flux( Ctx *, Vec );
PetscErrorCode set_mesh( Ctx * );
PetscErrorCode set_core_cooling( Ctx * );
PetscErrorCode set_twophase( Ctx * );

PetscScalar get_val1d( Interp1d *, PetscScalar );
PetscScalar get_val2d( Interp2d *, PetscScalar, PetscScalar );

PetscScalar combine_matprop( PetscScalar, PetscScalar, PetscScalar);

PetscScalar radiative_flux_with_dT( PetscScalar );
#endif
