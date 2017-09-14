#ifndef GLOBAL_DEFS_H_
#define GLOBAL_DEFS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// TODO : remove as much as possible from file, preferring to set things in parameters.c

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

/* number of additional equations for the augmented system */
//#define AUG_NUM 1 // for no coupled atmosphere evolution
// FIXME: this must always be set to 3
#define AUG_NUM 3 // for coupled atmosphere evolution

/* scaling constants */
static const PetscScalar RADIUS0 = 6371000.0; // m
static const PetscScalar ENTROPY0 = 2993.025100070677; // J/kg K
static const PetscScalar TEMPERATURE0 = 4033.6070755893948; // K
static const PetscScalar DENSITY0 = 4613.109568155063; // kg/m^3

/* datafile locations and material-specific constants */

/* A Root Directory */
/* https://gcc.gnu.org/onlinedocs/gcc-4.9.0/cpp/Stringification.html */
#define STRINGIFY(x) STRINGIFY2(x)
#define STRINGIFY2(x) #x
#define MAGMA_ROOT_DIR_STR STRINGIFY(MAGMA_ROOT_DIR)

/* liquidus/solidus data files */
#define LIQUIDUS_STIXRUDE2009 "../../../data/lookup/lookup-fusion/2009_stixrude/RTmelt/liquidus.dat"
#define LIQUIDUS_ANDRAULT2011 "../../../data/lookup/lookup-fusion/2011_andrault/RTmelt/liquidus.dat"
#define SOLIDUS_STIXRUDE2009 "../../../data/lookup/lookup-fusion/2009_stixrude/RTmelt/solidus.dat"
#define SOLIDUS_ANDRAULT2011 "../../../data/lookup/lookup-fusion/2011_andrault/RTmelt/solidus.dat"

/* solid data files */
#define ALPHA_SOL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/thermal_exp_solid.dat"

#define CP_SOL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/heat_capacity_solid.dat"

#define DTDPS_SOL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/adiabat_temp_grad_solid.dat"

#define RHO_SOL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/density_solid.dat"

#define TEMP_SOL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/temperature_solid.dat"

/* melt data files */
#define ALPHA_MEL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/thermal_exp_melt.dat"

#define CP_MEL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/heat_capacity_melt.dat"

#define DTDPS_MEL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/adiabat_temp_grad_melt.dat"

#define RHO_MEL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/density_melt.dat"

#define TEMP_MEL_DEFAULT "../../../data/lookup/lookup-hires-RTmelt/evo/temperature_melt.dat"

#endif
