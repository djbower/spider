#ifndef GLOBAL_DEFS_H_
#define GLOBAL_DEFS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* gsl library for 1D and 2D interpolation */
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/* PETSc (for data types) */
#include "petsc.h"

/* Define to have verbose and/or debugging output */
#define VERBOSE 1
#define DEBUGOUTPUT 1

#define HEAD 4 /* number of header lines in datafile */
#define NROWS 10100 /* number of coordinates in datafiles */
/* N.B., NROWS = NX * NY */
#define NX 101 /* no. of x coordinates in datafiles */
#define NY 100 /* no. of y coordinates in datafiles */

/* number of basic mesh points */
#define NUMPTS 100
/* number of staggered mesh points */
#define NUMPTSS NUMPTS-1

/* no. of coordinates in liq and sol datafiles */
#define NLS 301

/* inital entropy of adiabat */
static const int SINIT = 3000.0;
//static const int SINIT = 2500.0;
/* this initial condition spans all melt fractions and is useful
   for testing */
//static const int SINIT = 2100.0;

/* constants */
static const PetscScalar RADOUT = 6371000.0; // m
static const PetscScalar RADIN = 3504050.0; // m
static const PetscScalar RHOS = 4078.95095544; // kg/m^3
static const PetscScalar BETA = 1.1115348931E-7; // 1/m
static const PetscScalar GRAIN = 1.0E-3; // m
static const PetscScalar GRAVITY = -10.0; // m/s^2
static const PetscScalar F_THRESHOLD = 0.6; // non-dim
static const PetscScalar DF_TRANSITION = 0.15; // non-dim

/* for radiative thermal boundary condition at the top surface */
static const PetscScalar CONSTBC = 0.0027401185339112;
static const PetscScalar EXPBC = 1.6357699211550170;
static const PetscScalar SIGMA = 5.670373E-8; // Stefan-Boltzmann constant
static const PetscScalar TEQM = 273.0; // equilibrium temp of planet

/* for core-cooling */
static const PetscScalar CMBBC = 0.99;

/* end of constants */

/* datafile locations and material-specific constants */

/* liquidus data file */
static const char LIQUIDUS[] = "../../../data/lookup/lookup-hires-RTmelt/liquidus.dat";

/* solidus data file */
static const char SOLIDUS[] = "../../../data/lookup/lookup-hires-RTmelt/solidus.dat";

/* solid data files */
static const char ALPHA_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/thermal_exp_solid.dat";

static const char CP_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/heat_capacity_solid.dat";

static const char DTDPS_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/adiabat_temp_grad_solid.dat";

static const char RHO_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/density_solid.dat";

static const char TEMP_SOL[] = "../../../data/lookup/lookup-hires-RTmelt/temperature_solid.dat";

static const PetscScalar LOG10VISC_SOL = 21.0;

static const PetscScalar COND_SOL = 4.0;

/* melt data files */
static const char ALPHA_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/thermal_exp_melt.dat";

static const char CP_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/heat_capacity_melt.dat";

static const char DTDPS_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/adiabat_temp_grad_melt.dat";

static const char RHO_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/density_melt.dat";

static const char TEMP_MEL[] = "../../../data/lookup/lookup-hires-RTmelt/temperature_melt.dat";

static const PetscScalar LOG10VISC_MEL = 2.0;

static const PetscScalar COND_MEL = 4.0;

#endif
