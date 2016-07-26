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

/* set mesh here */
/* number of basic mesh points */
#define NUMPTS 100
/* number of staggered mesh points */
#define NUMPTSS NUMPTS-1 /* automagically determined */

/* initial condition: set entropy of adiabat */
//static const PetscInt SINIT = 3000.0; // all super-liquidus
//static const PetscInt SINIT = 1600.0; // all sub-solidus
/* this initial condition spans all melt fractions and is useful
   for testing to make sure the RHS is correct */
static const PetscInt SINIT = 2500.0;
/* percent linear increase in entropy with pressure to make the
   initial condition just slightly super-adiabatic */
static const PetscScalar SUPERFAC = 0.01;

/* constants */
/* outer radius */
static const PetscScalar RADOUT = 6371000.0; // m
/* inner radius */
static const PetscScalar RADIN = 3504050.0; // m
/* surface density for Adams-Williamson EOS for pressure */
static const PetscScalar RHOS = 4078.95095544; // kg/m^3
/* parameter for Adams-Williamson EOS for pressure */
static const PetscScalar BETA = 1.1115348931E-7; // 1/m
/* grain size */
static const PetscScalar GRAIN = 1.0E-3; // m
/* gravity.  Always constant and must be negative */
static const PetscScalar GRAVITY = -10.0; // m/s^2
/* melt fraction threshold for rheology */
static const PetscScalar F_THRESHOLD = 0.6; // non-dim
/* melt fraction transition width for rheology */
static const PetscScalar DF_TRANSITION = 0.15; // non-dim
/* melt fraction shape transition for skew */
static const PetscScalar SHAPE_TRANSITION = 0.5;

/* for radiative thermal boundary condition at the top surface */
/* dT = CONSTBC * [Potential temperature]**EXPBC */
static const PetscScalar CONSTBC = 0.0027401185339112;
static const PetscScalar EXPBC = 1.6357699211550170;
static const PetscScalar SIGMA = 5.670373E-8; // Stefan-Boltzmann constant
static const PetscScalar TEQM = 273.0; // equilibrium temp of planet

/* for core-cooling boundary condition at the bottom surface */
/* mass of core (perhaps excluding inner core?) */
static const PetscScalar MCORE = 1.9352467333195415E24; // kg
/* mass-weighted average core temperature as a fraction of CMB temp */
static const PetscScalar TFAC_CORE_AVG = 1.147;
/* heat capacity of core */
static const PetscScalar CP_CORE = 880.0; // # J/K/kg
/* density near CMB. TODO: this is f(S,P) and thus time-dep */
static const PetscScalar RHO_CMB = 5500.0; // kg/m^3
/* heat capacity near CMB.  TODO: this is f(S,P) and thus time-dep */
static const PetscScalar CP_CMB = 1200.0; // J/K/kg

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
