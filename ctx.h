#ifndef CTX_H_
#define CTX_H_

#include "petsc.h"

// Heinous. Petsc likes include headers which define I as a complex number,
// even if you aren't aware that PETSc even supports that. 
#undef I

#include "global_defs.h"

// !! Temp until we move the Ctx constructor and destructor to another file
#include "ctx.h"

/* common structures */

/* this structure has an x and y array size equal to NLS */
typedef struct _Interp1d {
    PetscScalar xa[NLS];
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar ya[NLS];
    PetscScalar ymin;
    PetscScalar ymax;
    gsl_interp *interp;
    gsl_interp_accel *acc;
} Interp1d;

typedef struct _Interp2d {
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar ymin;
    PetscScalar ymax;
    gsl_spline2d *interp;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
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

#define NUMMESHVECS 5
#define NUMMESHVECSS 3
typedef struct _Mesh {
    PetscScalar dx_b;
    PetscScalar dx_s;

    Vec meshVecs[NUMMESHVECS];
    Vec area_b,dPdr_b,pressure_b,radius_b, mix_b;  

    Vec meshVecsS[NUMMESHVECSS];
    Vec pressure_s,radius_s,volume_s;

} Mesh;

#define NUMSOLUTIONVECS 40
#define NUMSOLUTIONVECSS 16
typedef struct _Solution {
    /* TI means time-independent */

    Vec solutionVecs[NUMSOLUTIONVECS];
    Vec alpha,alpha_mix,cond,cp,cp_mix,dfusdr,dfusdr_temp,dphidr,dSdr,dTdPs,dTdPs_mix,dTdrs,dTdrs_mix,Etot,fusion,fusion_curve,fusion_curve_temp,fusion_rho,fusion_temp,gsuper,Jcond,Jconv,Jgrav,Jheat,Jmass,Jmix,Jtot,kappah,liquidus,liquidus_rho,liquidus_temp,nu,phi,rho,S,solidus,solidus_rho,solidus_temp,temp,visc;

    Vec solutionVecsS[NUMSOLUTIONVECSS];
    Vec fusion_s, fusion_curve_s, fusion_curve_temp_s, fusion_temp_s, lhs_s, liquidus_rho_s, liquidus_s, liquidus_temp_s, phi_s, rhs_s, rho_s, S_s, solidus_s, solidus_rho_s, solidus_temp_s, temp_s; 

} Solution;

/* A Context for the Solver */
typedef struct _Ctx {
  // Here we essentially copy what is in the All_variables field
  // The items in the mesh and solution sub-structs there need to be replaced with versions
  //   based on vectors (we don't worry about AoS vs SoA at this point)
  //..
  Lookup   melt_prop;
  Lookup   solid_prop;
  Mesh     mesh;
  Solution solution;
  DM       da_b,da_s;
  PetscScalar BC_BOT_FAC; // core-cooling boundary condition

  /* "local" work vectors */
  Vec work_local_s,work_local_b;

  // TODO: flatten this all out
} Ctx;

PetscErrorCode d_dr( Ctx *, Vec, Vec );
PetscErrorCode free_memory_interp( Ctx * );
PetscErrorCode set_capacitance( Ctx *, Vec );
PetscErrorCode set_initial_condition( Ctx * );
PetscErrorCode set_lookups( Ctx * );
PetscErrorCode set_matprop_and_flux( Ctx * );
PetscErrorCode set_mesh( Ctx * );
PetscErrorCode set_time_independent( Ctx * );

PetscErrorCode set_interp2d( const char * filename, Interp2d *interp );
PetscErrorCode set_interp1d( const char * filename, Interp1d *interp, int n );

PetscScalar combine_matprop( PetscScalar, PetscScalar, PetscScalar);

PetscScalar average( PetscScalar a, PetscScalar b );

#endif
