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
#define GREYBODY
//#define HAMANO
//#define ZAHNLE

/* if HYBRID is set, then the boundary condition will switch
   to the upper mantle cooling rate once the rheological
   transition is reached.  This prevents a lid from forming at
   the top of the model.

   TODO: this implies that the emissivity is around 1.0E-7,
   which is unphysical unless you are appealing to a massive
   massive atmosphere */

#define HYBRID

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

/* number of additional equations for the augmented system */
//#define AUG_NUM 1 // for no coupled atmosphere evolution
// FIXME: this must always be set to 3
#define AUG_NUM 3 // for coupled atmosphere evolution

/* number of staggered mesh points */
#define NUMPTS_S_DEFAULT NUMPTS_B_DEFAULT-1 /* automagically determined */

/* dimensional constants */
/* FIXME: these constants must have been used to non-dimensionalise
   the data tables that are read in, since the header line is ignored
   for non-dimensional runs */
// must turn on LOOKUPDIM if you are using dimensional quantities
//#define LOOKUPDIM 1
static const PetscScalar RADIUS0 = 6371000.0; // m
static const PetscScalar ENTROPY0 = 2993.025100070677; // J/kg K
static const PetscScalar TEMPERATURE0 = 4033.6070755893948; // K
static const PetscScalar DENSITY0 = 4613.109568155063; // kg/m^3

/* initial condition: set entropy of adiabat */
/* set at reference entropy */
/* to compare with python test1 */
//static const PetscScalar SINIT_DEFAULT = 3052.885602072091;
static const PetscScalar SINIT_DEFAULT = 1.02;

/* to compare with python test2 */
/* TODO: refresh these tests */
//static const PetscScalar SINIT_DEFAULT = 0.7;

/* to compare with python test3 */
//static const PetscScalar SINIT_DEFAULT = 0.5;

/* initial entropy gradient */
//static const PetscScalar IC_DSDR = -4.6978890285209187e-07;
static const PetscScalar IC_DSDR = -1.0E-3;

/* planetary radius */
//static const PetscScalar RADIUS = 6371000.0; // m (outer radius of planet)
//static const PetscScalar RADIUS = 3185500.0; // m (outer radius of planet)
static const PetscScalar RADIUS = 1.0;

/* seconds in 1 year */
/* FIXME: scaling the rhs by this value could blow up noise */
//static const PetscScalar SECSINYR = 31557600.0; // s

/* core-mantle boundary radius */
static const PetscScalar CORESIZE = 0.55; // radius fraction

/* surface density for Adams-Williamson EOS for pressure */
//static const PetscScalar RHOS = 4078.95095544;
static const PetscScalar RHOS = 0.8842085571948184;
/* parameter for Adams-Williamson EOS for pressure */
//static const PetscScalar BETA = 1.1115348931000002e-07; // m^{-1}
static const PetscScalar BETA = 0.7081588803940101;
/* grain size */
//static const PetscScalar GRAIN = 1.0e-3;
static const PetscScalar GRAIN = 1.56961230576e-10;
/* gravity.  Always constant and must be negative */
//static const PetscScalar GRAVITY = -10.0;
static const PetscScalar GRAVITY = -5.2772012422265826;
/* melt fraction threshold for rheology */
static const PetscScalar PHI_CRITICAL = 0.4;
/* melt fraction transition width for rheology */
/* when PHI_WIDTH = 0.15, oscillations appear in the volatile contents
   due to a rigid crust forming at the top of the model.  Reducing the value
   to 0.2 helps to alleviate this problem.  So evidently the viscosity contrast
   across nodes matters. */
static const PetscScalar PHI_WIDTH = 0.15;
// testing atmosphere evolution and boundary layer formation
//static const PetscScalar PHI_WIDTH = 0.2;
/* melt fraction shape transition for skew */
static const PetscScalar PHI_SKEW = 0.0;

/* for radiative thermal boundary condition at the top surface */
/* dT = CONSTBC * [Surface temperature]**3 */
/* FIXME: I think in this form you have to rescale by TEMPERATURE0
   in the code.  To check */
//static const PetscScalar CONSTBC = 1.0e-07; // 1.0E-7 in python
static const PetscScalar CONSTBC = 1.6269986040244828; // 1.0E-7 in python
//static const PetscScalar CONSTBC = 0.0; // no ultra-thin thermal boundary layer

/* for core-cooling boundary condition at the bottom surface */
/* density of core */
//static const PetscScalar RHO_CORE = 10738.332568062382;
static const PetscScalar RHO_CORE = 2.3277861514910865;
/* heat capacity of core */
//static const PetscScalar CP_CORE = 880.0;
static const PetscScalar CP_CORE = 0.2940169128482149;
/* mass-weighted average core temperature as a fraction of CMB temp */
static const PetscScalar TFAC_CORE_AVG = 1.147;

/* smoothing parameters */
/* width */
static const PetscScalar SWIDTH = 1.0E-2;

/* atmosphere parameters */
/* NOTE: if both H2O_INITIAL and CO2_INITIAL are set to zero, then
   the emissivity is constant with time (a grey-body, i.e.,
   a black-body scaled by EMISSIVITY.  If either H2O_INITIAL and/or
   CO2_INITIAL are positive then the emissivity is computed according
   to the plane-parallel radiative equilibrium model of Abe and 
   Matsui (1985)  */
// volatile scaling
/* VOLSCALE enables us to scale the volatile equations to the same
   order of magnitude as entropy, and thus ensure that the residual
   based on the solution vector is not biased. */
//static const PetscScalar VOLSCALE = 100.0 // for wt %.
static const PetscScalar VOLSCALE = 1000000.0; // for ppm

// Elkins_Tanton (2008) estimates of volatile content
//static const PetscScalar H2O_INITIAL = 0.0; // units according to VOLSCALE
// Turn off atmosphere
static const PetscScalar H2O_INITIAL = 0.0;
static const PetscScalar CO2_INITIAL = 0.0;
// ET08 Case 1
//static const PetscScalar H2O_INITIAL = 500.0;
//static const PetscScalar CO2_INITIAL = 100.0;
// ET08 Case 2
//static const PetscScalar H2O_INITIAL = 5000.0;
//static const PetscScalar CO2_INITIAL = 1000.0;
// ET08 Case 3
//static const PetscScalar H2O_INITIAL = 0.0;
//static const PetscScalar CO2_INITIAL = 6000.0;

// EMISSIVITY is not used if H20_INITIAL and/or CO2_INITIAL > 0.0
static const PetscScalar EMISSIVITY = 1.0;

//static const PetscScalar SIGMA = 5.670367e-08; // Stefan-Boltzmann constant
static const PetscScalar SIGMA = 7.75685789946723e-08;

//static const PetscScalar TEQM = 273.0; // equilibrium temp of planet
static const PetscScalar TEQM = 0.0676813568808283;

#if 1
// Elkins-Tanton (2008) parameters
static const PetscScalar P0 = 101325.0; // Pa (= 1 atm)
// distribution coefficients are given in supplementary information
static const PetscScalar H2O_KDIST = 1.0E-4; // distribution coefficient between solid and melt
// TODO: water saturation limit is 10 ppm
static const PetscScalar H2O_KABS = 0.01; // m^2/kg
// next two from Lebrun et al. (2013)
static const PetscScalar H2O_HENRY = 6.8E-8; // must be mass fraction/Pa
static const PetscScalar H2O_HENRY_POW = 1.4285714285714286; // (1.0/0.7)

static const PetscScalar CO2_KDIST = 5.0E-4; // distribution coefficient between solid and melt
// TODO: carbon saturation limit is 0.03 ppm
static const PetscScalar CO2_KABS = 0.05; // m^2/kg
//static const PetscScalar CO2_KABS = 0.0001; // m^2/kg
// CO2 dissolved in silicate melt obeys Henry's law to a good approximation
// e.g. Pan et al., 1991
// CO2_HENRY = concentration (mass fraction) / partial pressure (Pa)
// next two from Lebrun et al. (2013)
static const PetscScalar CO2_HENRY = 4.4E-12; // must be mass fraction/Pa
static const PetscScalar CO2_HENRY_POW = 1.0;
#endif

#if 0
// Salvador et al. (2017)
// TODO: check these parameters
static const PetscScalar P0 = 101325.0; // Pa

static const PetscScalar H2O_KDIST = 1.0E-4; // (from ET08)
static const PetscScalar H2O_INITIAL = 1.17E-3; // wt. %
//static const PetscScalar H2O_INITIAL = 4.67E-2;
static const PetscScalar H2O_KABS = 0.01; // m^2/kg

static const PetscScalar CO2_KDIST = 5.0E-4; // (from ET08)
static const PetscScalar CO2_INITIAL = 1.0E-3;
static const PetscScalar CO2_INITIAL = 1.4E-1;
static const PetscScalar CO2_KABS = 0.0001; // m^2/kg
#endif

/* end of constants */

/* datafile locations and material-specific constants */

/* liquidus data file */
static const char LIQUIDUS[] = "../../../data/lookup/lookup-fusion/2009_stixrude/RTmelt/liquidus.dat";
//static const char LIQUIDUS[] = "../../../data/lookup/lookup-fusion/2011_andrault/RTmelt/liquidus.dat";

/* solidus data file */
static const char SOLIDUS[] = "../../../data/lookup/lookup-fusion/2009_stixrude/RTmelt/solidus.dat";
//static const char SOLIDUS[] = "../../../data/lookup/lookup-fusion/2011_andrault/RTmelt/solidus.dat";

/* solid data files */
static const char ALPHA_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/thermal_exp_solid.dat";

static const char CP_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/heat_capacity_solid.dat";

static const char DTDPS_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/adiabat_temp_grad_solid.dat";

static const char RHO_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/density_solid.dat";

static const char TEMP_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/temperature_solid.dat";

//static const PetscScalar LOG10VISC_SOL = 21.0;
static const PetscScalar LOG10VISC_SOL = 6.99089665051;

//static const PetscScalar COND_SOL = 4.0;
static const PetscScalar COND_SOL = 1.30871862439e-17;

/* melt data files */
static const char ALPHA_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/thermal_exp_melt.dat";

static const char CP_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/heat_capacity_melt.dat";

static const char DTDPS_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/adiabat_temp_grad_melt.dat";

static const char RHO_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/density_melt.dat";

static const char TEMP_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/evo/temperature_melt.dat";

//static const PetscScalar LOG10VISC_MEL = 2.0;
static const PetscScalar LOG10VISC_MEL = -12.0091033495;

//static const PetscScalar COND_MEL = 4.0;
static const PetscScalar COND_MEL = 1.30871862439e-17;

#endif
