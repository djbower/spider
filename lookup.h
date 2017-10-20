#ifndef LOOKUP_H_
#define LOOKUP_H_

#include <petsc.h>

/* 2-D datafiles containing melt and solid properties
   as a functon of pressure and entropy */
#define HEAD 4 /* number of header lines in datafile */

/* these next values are for lookup-hires-{RT}melt */
/* note that you must also change the file locations below! */
#define NROWS 10100 /* number of coordinates in datafiles */
/* N.B., NROWS = NX * NY */
#define NX 101 /* no. of x coordinates in datafiles */
#define NY 100 /* no. of y coordinates in datafiles */

/* these next values are for lookup-hires-{RT}press */
//#define NROWS 15100 /* number of coordinates in datafiles */
/* N.B., NROWS = NX * NY */
//#define NX 151 /* no. of x coordinates in datafiles */
//#define NY 100 /* no. of y coordinates in datafiles */

/* 1-D datafiles containing liquidus and solidus
   as a function of pressure */
#define NLS 301 /* no. of coordinates in liq and sol datafiles */

/* datafile locations and material-specific constants */

/* A Root Directory */
/* https://gcc.gnu.org/onlinedocs/gcc-4.9.0/cpp/Stringification.html */
#define STRINGIFY(x) STRINGIFY2(x)
#define STRINGIFY2(x) #x
#define MAGMA_ROOT_DIR_STR STRINGIFY(MAGMA_ROOT_DIR)

/* liquidus/solidus data files */
#define LIQUIDUS_STIXRUDE2009 "lookup_data/RTmelt/liquidus_stixrude2009.dat"
#define LIQUIDUS_ANDRAULT2011 "lookup_data/RTmelt/liquidus_andrault2011.dat"
#define SOLIDUS_STIXRUDE2009 "lookup_data/RTmelt/solidus_stixrude2009.dat"
#define SOLIDUS_ANDRAULT2011 "lookup_data/RTmelt/solidus_andrault2011.dat"

/* solid data files */
#define ALPHA_SOL_DEFAULT "lookup_data/RTmelt/thermal_exp_solid.dat"
#define CP_SOL_DEFAULT "lookup_data/RTmelt/heat_capacity_solid.dat"
#define DTDPS_SOL_DEFAULT "lookup_data/RTmelt/adiabat_temp_grad_solid.dat"
#define RHO_SOL_DEFAULT "lookup_data/RTmelt/density_solid.dat"
#define TEMP_SOL_DEFAULT "lookup_data/RTmelt/temperature_solid.dat"

/* melt data files */
#define ALPHA_MEL_DEFAULT "lookup_data/RTmelt/thermal_exp_melt.dat"
#define CP_MEL_DEFAULT "lookup_data/RTmelt/heat_capacity_melt.dat"
#define DTDPS_MEL_DEFAULT "lookup_data/RTmelt/adiabat_temp_grad_melt.dat"
#define RHO_MEL_DEFAULT "lookup_data/RTmelt/density_melt.dat"
#define TEMP_MEL_DEFAULT "lookup_data/RTmelt/temperature_melt.dat"

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
