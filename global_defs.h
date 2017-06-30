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

/* select which energy terms to include */
#define CONDUCTION
#define CONVECTION
#define MIXING
#define SEPARATION
//#define HRADIO
//#define HTIDAL

/* select which surface boundary condition to use */
/* these are mutually exclusive, so only leave one
   of them uncommented */
//#define GREYBODY
//#define HAMANO
//#define HYBRID
#define ZAHNLE

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
///static const PetscScalar SINIT_DEFAULT = 1.02;
static const PetscScalar SINIT_DEFAULT = 3052.885602072091;

/* to compare with python test2 */
/* TODO: refresh these tests */
//static const PetscScalar SINIT_DEFAULT = 0.7;

/* to compare with python test3 */
//static const PetscScalar SINIT_DEFAULT = 0.5;

/* initial entropy gradient */
///static const PetscScalar IC_DSDR = -1.0E-3;
static const PetscScalar IC_DSDR = -4.6978890285209187e-07;

/* constants */
/* do not change any of the next four constants without first
   consulting Dan! */
static const PetscScalar RADOUT = 6371000.0; // m (outer radius of planet)
///static const PetscScalar ENTROPY0 = 2993.025100070677; // J/kg K
///static const PetscScalar TEMPERATURE0 = 4033.6070755893948; // K
///static const PetscScalar DENSITY0 = 4613.109568155063; // kg/m^3

/* DJB TO SORT OUT DIMENSIONAL SCALINGS BELOW */
/* dimensional time scaling */
///#define TIME0YEARS RADIUS0/PetscSqrtScalar(ENTROPY0*TEMPERATURE0)/(3600.0*24.0*365.25); // 5.810341565106721e-05;

/* seconds in 1 year */
static const PetscScalar SECSINYR = 31557600.0; // s

/* dimensional flux scaling */
///#define FLUX0 DENSITY0*PetscPowScalar(ENTROPY0*TEMPERATURE0,3.0/2.0);

/* core-mantle boundary radius */
///static const PetscScalar RADIN = 0.55;
static const PetscScalar RADIN = 3504050.0; // m (inner radius of mantle)
/* surface density for Adams-Williamson EOS for pressure */
//static const PetscScalar RHOS = 0.8842085571948184;
static const PetscScalar RHOS = 4078.95095544;
/* parameter for Adams-Williamson EOS for pressure */
//static const PetscScalar BETA = 0.7081588803940101;
static const PetscScalar BETA = 1.1115348931000002e-07; // m^{-1}
/* grain size */
static const PetscScalar GRAIN = 1.0e-3;
/* gravity.  Always constant and must be negative */
static const PetscScalar GRAVITY = -10.0;
/* melt fraction threshold for rheology */
static const PetscScalar PHI_CRITICAL = 0.4;
/* melt fraction transition width for rheology */
static const PetscScalar PHI_WIDTH = 0.15;
/* melt fraction shape transition for skew */
static const PetscScalar PHI_SKEW = 0.0;

/* for radiative thermal boundary condition at the top surface */
/* dT = CONSTBC * [Surface temperature]**3 */
/* where dT and Surface temperature are non-dimensional */
/* you can extract the value of CONSTBC by looking at the non-dimensional
   value in a run.log output by the python evolution code.  Or alternatively,
   CONSTBC here is [dimensional CONSTBC] * [TEMP0]**2, where TEMP0 is the
   dimensional (scaling) temperature.  See python script for more information */
//static const PetscScalar CONSTBC = 0.0;
//static const PetscScalar CONSTBC = 1.6269986040244828; // 1.0E-7 in python
static const PetscScalar CONSTBC = 1.0e-07; // 1.0E-7 in python
/* SB constant non-dimensionalisation depends on RADIUS0 */
//static const PetscScalar SIGMA = 7.75685789946723e-08; // Stefan-Boltzmann constant
static const PetscScalar SIGMA = 5.670367e-08; // Stefan-Boltzmann constant
static const PetscScalar EMISSIVITY = 1.0; // emissivity
///static const PetscScalar TEQM = 0.0676813568808283; // equilibrium temp of planet
static const PetscScalar TEQM = 273.0; // equilibrium temp of planet

/* for core-cooling boundary condition at the bottom surface */
/* density of core */
///static const PetscScalar RHO_CORE = 2.3277861514910865;
static const PetscScalar RHO_CORE = 10738.332568062382;
/* heat capacity of core */
///static const PetscScalar CP_CORE = 0.2940169128482149;
static const PetscScalar CP_CORE = 880.0;
/* mass-weighted average core temperature as a fraction of CMB temp */
static const PetscScalar TFAC_CORE_AVG = 1.147;

/* smoothing parameters */
/* width */
///static const PetscScalar SWIDTH = 1.0E-2;
/* these units stay the same (units of melt fraction) */
static const PetscScalar SWIDTH = 1.0E-2;

/* end of constants */

/* datafile locations and material-specific constants */

/* liquidus data file */
//static const char LIQUIDUS[] = "../../../data/lookup/lookup-fusion/2009_stixrude/RTmelt/liquidus.dat";
static const char LIQUIDUS[] = "../../../data/lookup/lookup-fusion/2011_andrault/RTmelt/liquidus.dat";

/* solidus data file */
//static const char SOLIDUS[] = "../../../data/lookup/lookup-fusion/2009_stixrude/RTmelt/solidus.dat";
static const char SOLIDUS[] = "../../../data/lookup/lookup-fusion/2011_andrault/RTmelt/solidus.dat";

/* solid data files */
static const char ALPHA_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/thermal_exp_solid.dat";

static const char CP_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/heat_capacity_solid.dat";

static const char DTDPS_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/adiabat_temp_grad_solid.dat";

static const char RHO_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/density_solid.dat";

static const char TEMP_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/temperature_solid.dat";

static const PetscScalar LOG10VISC_SOL = 21.0;

static const PetscScalar COND_SOL = 4.0;

/* melt data files */
static const char ALPHA_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/thermal_exp_melt.dat";

static const char CP_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/heat_capacity_melt.dat";

static const char DTDPS_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/adiabat_temp_grad_melt.dat";

static const char RHO_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/density_melt.dat";

static const char TEMP_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/temperature_melt.dat";

static const PetscScalar LOG10VISC_MEL = 2.0;

static const PetscScalar COND_MEL = 4.0;

#endif
