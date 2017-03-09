#ifndef GLOBAL_DEFS_H_
#define GLOBAL_DEFS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* PETSc (for data types) */
#include "petsc.h"

/* Definitions */
//#define VERBOSE 1
//#define DEBUGOUTPUT 1

/* 2-D datafiles containing melt and solid properties
   as a functon of pressure and entropy */
#define HEAD 4 /* number of header lines in datafile */
#define NROWS 10100 /* number of coordinates in datafiles */
/* N.B., NROWS = NX * NY */
#define NX 101 /* no. of x coordinates in datafiles */
#define NY 100 /* no. of y coordinates in datafiles */

/* 1-D datafiles containing liquidus and solidus
   as a function of pressure */
#define NLS 301 /* no. of coordinates in liq and sol datafiles */

/* set default mesh here (can be changed from the command line) */
/* number of basic mesh points */
/* you have to use the python script to a priori determine how many
   nodes there are in a given mesh */

/* for constant mesh spacing */
#define NUMPTS_B_DEFAULT 200

//#define NUMPTS_B_DEFAULT 278
//#define NUMPTS_B_DEFAULT 372
//#define NUMPTS_B_DEFAULT 656
//#define NUMPTS_B_DEFAULT 2939
//#define NUMPTS_B_DEFAULT 5802
//#define NUMPTS_B_DEFAULT 11532
//#define NUMPTS_B_DEFAULT 22996
//#define NUMPTS_B_DEFAULT 45928

/* number of staggered mesh points */
#define NUMPTS_S_DEFAULT NUMPTS_B_DEFAULT-1 /* automagically determined */

/* initial condition: set entropy of adiabat */
/* set at reference entropy */
/* to compare with python test1 */
static const PetscScalar SINIT_DEFAULT = 1.02;

/* to compare with python test2 */
/* TODO: refresh these tests */
//static const PetscScalar SINIT_DEFAULT = 0.8;

/* to compare with python test3 */
//static const PetscScalar SINIT_DEFAULT = 0.5;

/* initial entropy gradient */
static const PetscScalar IC_DSDR = -1.0E-10;

/* constants */
/* core-mantle boundary radius */
static const PetscScalar RADIN = 0.55;
/* surface density for Adams-Williamson EOS for pressure */
static const PetscScalar RHOS = 0.8842085571948184;
/* parameter for Adams-Williamson EOS for pressure */
static const PetscScalar BETA = 0.7081588803940101;
/* grain size */
static const PetscScalar GRAIN = 1.56961230576e-10;
/* gravity.  Always constant and must be positive */
static const PetscScalar GRAVITY = 5.2772012422265826;
/* melt fraction threshold for rheology */
static const PetscScalar F_THRESHOLD = 0.6;
/* melt fraction transition width for rheology */
static const PetscScalar DF_TRANSITION = 0.15;
/* melt fraction shape transition for skew */
static const PetscScalar SHAPE_TRANSITION = 0.5;

/* for radiative thermal boundary condition at the top surface */
/* dT = CONSTBC * [Potential temperature]**EXPBC */
static const PetscScalar CONSTBC = 0.0; // 0.0036006849671025;
static const PetscScalar EXPBC = 0.0; //1.6042079063944077;
static const PetscScalar SIGMA = 7.75685789946723e-08; // Stefan-Boltzmann constant
static const PetscScalar EMISSIVITY = 1.0; // emissivity
static const PetscScalar TEQM = 0.0676813568808283; // equilibrium temp of planet

/* for core-cooling boundary condition at the bottom surface */
/* mass of core (perhaps excluding inner core?) */
static const PetscScalar MCORE = 1.62225737775007;
/* heat capacity of core */
static const PetscScalar CP_CORE = 0.2940169128482149;
/* mass-weighted average core temperature as a fraction of CMB temp */
static const PetscScalar TFAC_CORE_AVG = 1.147;

/* smoothing parameters */
/* width */
static const PetscScalar SWIDTH = 1.0E-2;

/* end of constants */

/* datafile locations and material-specific constants */

/* liquidus data file */
static const char LIQUIDUS[] = "../../../data/lookup/lookup-fusion/bridgmanite/liquidus.dat";

/* solidus data file */
static const char SOLIDUS[] = "../../../data/lookup/lookup-fusion/bridgmanite/solidus.dat";

/* solid data files */
static const char ALPHA_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/thermal_exp_solid.dat";

static const char CP_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/heat_capacity_solid.dat";

static const char DTDPS_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/adiabat_temp_grad_solid.dat";

static const char RHO_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/density_solid.dat";

static const char TEMP_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/temperature_solid.dat";

static const PetscScalar LOG10VISC_SOL = 6.99089665051;

static const PetscScalar COND_SOL = 1.30871862439e-17;

/* melt data files */
static const char ALPHA_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/thermal_exp_melt.dat";

static const char CP_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/heat_capacity_melt.dat";

static const char DTDPS_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/adiabat_temp_grad_melt.dat";

static const char RHO_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/density_melt.dat";

static const char TEMP_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/temperature_melt.dat";

static const PetscScalar LOG10VISC_MEL = -12.0091033495;

static const PetscScalar COND_MEL = 1.30871862439e-17;

#endif
