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
    PetscInt    n;
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

/* for storing atmosphere outputs for eventual writing to Petsc
   binary file */
typedef enum {MO_ATMOSPHERE_TYPE_GREY_BODY=1,MO_ATMOSPHERE_TYPE_ZAHNLE,MO_ATMOSPHERE_TYPE_VOLATILES} MagmaOceanAtmosphereType;
typedef struct _Atmosphere {
    // 37 parameters in total
    // input parameters (20)
    MagmaOceanAtmosphereType MODEL;
    PetscInt HYBRID;
    // below are standard, also used for grey-body atmosphere
    PetscScalar EMISSIVITY;
    PetscScalar SIGMA;
    PetscScalar TEQM;
    PetscScalar CONSTBC;
    // for volatile ODE
    PetscScalar VOLSCALE;
    PetscScalar P0;
    // H2O (TODO: move to a 'Volatile' struct)
    PetscScalar H2O_INITIAL;
    PetscScalar H2O_KDIST;
    PetscScalar H2O_KABS;
    PetscScalar H2O_HENRY;
    PetscScalar H2O_HENRY_POW;
    // CO2 (TODO: move to a 'Volatile' struct)
    PetscScalar CO2_INITIAL;
    PetscScalar CO2_KDIST;
    PetscScalar CO2_KABS;
    PetscScalar CO2_HENRY;
    PetscScalar CO2_HENRY_POW;
    /* although RADIUS and GRAVITY are duplicated, they may have
       different non-dimensional values depending on the volatile
       non-dimenisonalisation scheme, which will be different to
       the dS/dr scheme to ensure comparable magnitude residuals */
    PetscScalar RADIUS; // duplicate
    PetscScalar GRAVITY; // duplicate
    // calculated quantites (17)
    PetscScalar M0; // total mass of mantle from EOS (kg)
    PetscScalar Mliq; // mass of liquid (kg)
    PetscScalar Msol; // mass of solid (kg)
    PetscScalar dMliqdt; // dMliq/dt (kg/yr)
    PetscScalar tau; // aggregate optical depth (dimensionless)
    PetscScalar x0; // CO2 content (wt %)
    PetscScalar dx0dt; // dx0/dt (wt % / yr)
    PetscScalar p0; // CO2 partial pressure (Pa)
    PetscScalar dp0dx; // dp0/dx (Pa/mass fraction)
    PetscScalar m0; // CO2 mass in atmosphere (kg)
    PetscScalar tau0; // CO2 optical depth (dimensionless)
    PetscScalar x1; // H2O content (wt %)
    PetscScalar dx1dt; // dx1/dt (wt % / yr)
    PetscScalar p1; // H2O partial pressure (Pa)
    PetscScalar dp1dx; // dp1dx (Pa / mass fraction)
    PetscScalar m1; // H2O mass in atmosphere (kg)
    PetscScalar tau1; // H20 optical depth (dimensionless)
} Atmosphere;

/* dimensionalising constants */
typedef struct _Constants {
    PetscScalar RADIUS;
    PetscScalar TEMP;
    PetscScalar ENTROPY;
    PetscScalar DENSITY;
    PetscScalar AREA;
    PetscScalar AREAG; // with 4*pi geometry
    PetscScalar VOLUME;
    PetscScalar VOLUMEG; // with 4*pi geometry
    PetscScalar MASS;
    PetscScalar MASSG; // with 4*pi geometry
    PetscScalar TIME;
    PetscScalar TIMEYRS;
    PetscScalar SENERGY;
    PetscScalar ENERGY;
    PetscScalar ENERGYG; // with 4*pi geometry
    PetscScalar PRESSURE;
    PetscScalar POWER;
    PetscScalar POWERG; // with 4*pi geometry
    PetscScalar FLUX;
    PetscScalar DPDR;
    PetscScalar ALPHA;
    PetscScalar GRAVITY;
    PetscScalar KAPPA;
    PetscScalar DTDP;
    PetscScalar DSDR;
    PetscScalar DTDR;
    PetscScalar GSUPER;
    PetscScalar ETA;
    PetscScalar LOG10ETA;
    PetscScalar NU;
    PetscScalar COND;
    PetscScalar SIGMA;
    PetscScalar LHS;
    PetscScalar LHSG; // with 4*pi geometry
    PetscScalar RHS;
} Constants;

typedef struct _Parameters {
    // 19
    PetscScalar sinit;
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
    PetscScalar swidth;
    PetscScalar log10visc_sol;
    PetscScalar cond_sol;
    PetscScalar log10visc_mel;
    PetscScalar cond_mel;
} Parameters;


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
#define NUMMESHVECS_S 7
typedef struct _Mesh {

    Vec meshVecs_b[NUMMESHVECS_B];
    Vec area_b,dPdr_b,pressure_b,radius_b, mix_b;  

    Vec meshVecs_s[NUMMESHVECS_S];
    Vec pressure_s,radius_s,volume_s,dPdr_s,area_s,rho_s,mass_s;

    /* DJB atmosphere.  For seeing what the 'pressure' estimate of the
       mass is */
    PetscScalar mass0;

} Mesh;

#define NUMSOLUTIONVECS_B 40
#define NUMSOLUTIONVECS_S 24
typedef struct _Solution {

    Vec solutionVecs_b[NUMSOLUTIONVECS_B];
    Vec alpha, alpha_mix, cond, cp, cp_mix, dfusdr, dfusdr_temp, dSdr, dSliqdr, dSsoldr, dTdrs, dTdrs_mix, Etot, fusion, fusion_curve, fusion_curve_temp, fusion_rho, fusion_temp, fwtl, fwts, gphi, gsuper, Jcond, Jconv, Jgrav, Jmix, Jtot, kappah, liquidus, liquidus_rho, liquidus_temp, nu, phi, rho, S, solidus, solidus_rho, solidus_temp, temp, visc;

    Vec solutionVecs_s[NUMSOLUTIONVECS_S];
    Vec cp_s, cp_mix_s, dSdt_s, fusion_s, fusion_curve_s, fusion_curve_temp_s, fusion_temp_s, fwtl_s, fwts_s, gphi_s, Hradio_s, Htidal_s, Htot_s, lhs_s, liquidus_rho_s, liquidus_s, liquidus_temp_s, phi_s, rho_s, S_s, solidus_s, solidus_rho_s, solidus_temp_s, temp_s;

} Solution;

/* A Context for the Solver */
typedef struct _Ctx {
  Lookup   melt_prop;
  Lookup   solid_prop;
  Mesh     mesh;
  Solution solution;
  DM       da_b,da_s;
  PetscScalar S_init; // initial entropy
  Mat      qty_at_b, ddr_at_b;
  Atmosphere atmosphere;
  Constants constants;
  Parameters parameters;

  /* "local" work vectors */
  Vec work_local_s,work_local_b;

  // Note: there is no real reason that we have some data here and some in Solution and some in Mesh.
  //       it could all be flattened out (and probably would be in a re-implementation)
} Ctx;

/* prototype definitions must come at the end, since before this
   point the code doesn't know about the Ctx structure */
PetscErrorCode setup_ctx(Ctx *);
PetscErrorCode destroy_ctx(Ctx *);

#endif
