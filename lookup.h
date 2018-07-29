#ifndef LOOKUP_H_
#define LOOKUP_H_

#include <petsc.h>

/* 2-D datafiles containing melt and solid properties
   as a functon of pressure and entropy */
#define HEAD 4 /* number of header lines in datafiles */

#if 0
/*------------------- */
/* lookup_data/RTmelt */
/* number of coordinates in datafiles */
#define NROWS 10100 /* RTmelt */
/* N.B., NROWS = NX * NY */
#define NX 101 /* no. of x coordinates in datafiles */
#define NY 100 /* no. of y coordinates in datafiles */
/* 1-D datafiles containing liquidus and solidus
   as a function of pressure */
#define NLS 301 /* no. of coordinates in liq and sol datafiles */
/*--------------------*/
#endif

#if 1
/*-------------------------------- */
/* lookup_data/1TPa-dK09-elec-free */
/* FIXME: these are actually for solid data tables only */
/* FIXME: melt tables have different length */
/* number of coordinates in datafiles */
#define NROWS 230280 /* solid only */ 
/* N.B., NROWS = NX * NY */
#define NX 2020 /* no. of x coordinates in datafiles */
#define NY 114 /* no. of y coordinates in datafiles */
/* 1-D datafiles containing liquidus and solidus
   as a function of pressure */
#define NLS 301 /* no. of coordinates in liq and sol datafiles */
/*---------------------------------*/
#endif


/* this structure has an x and y array size equal to NLS  */
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

PetscErrorCode set_interp2d( const char *, Interp2d *, PetscScalar, PetscScalar, PetscScalar );
PetscErrorCode set_interp1d( const char *, Interp1d *, PetscInt, PetscScalar, PetscScalar );
PetscScalar get_val1d( Interp1d const *, PetscScalar );
PetscScalar get_val2d( Interp2d const *, PetscScalar, PetscScalar );

#endif
