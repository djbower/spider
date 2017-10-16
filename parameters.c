/*
Parameter Management

Parameters should only ever be set by the functions in this file. That is, everywhere else they should be considered read-only.

Custom PETSc command line options should only ever be parsed here.
*/

#include "parameters.h"
#include "ctx.h"

static PetscErrorCode SetConstants( Constants *C, PetscReal RADIUS, PetscReal TEMPERATURE, PetscReal ENTROPY, PetscReal DENSITY )
{
    PetscScalar SQRTST;

    PetscFunctionBeginUser;
    /* 28 constants to set (excluding SQRTST which is a convenience
       parameter) */
    SQRTST = PetscSqrtScalar( ENTROPY * TEMPERATURE );

    C->RADIUS    = RADIUS; // m
    C->TEMP      = TEMPERATURE; // K
    C->ENTROPY   = ENTROPY; // (specific) J/kg.K
    C->DENSITY   = DENSITY; // kg/m^3
    C->AREA      = PetscSqr( C->RADIUS ); // m^2
    C->VOLUME    = C->AREA * C->RADIUS; // m^3
    C->MASS      = C->DENSITY * C->VOLUME; // kg
    C->TIME      = C->RADIUS / SQRTST; // s
    C->TIMEYRS   = C->TIME / (60.0*60.0*24.0*365.25); // years
    C->SENERGY   = C->ENTROPY * C->TEMP; // J/kg
    C->ENERGY    = C->SENERGY * C->MASS; // J
    C->PRESSURE  = C->ENTROPY * C->TEMP * C->DENSITY; // Pa
    C->POWER     = C->ENERGY / C->TIME; // W
    C->FLUX      = C->POWER / C->AREA; // W/m^2
    C->DPDR      = C->PRESSURE / C->RADIUS; // Pa/m
    C->GRAVITY   = (C->ENTROPY * C->TEMP) / C->RADIUS; // m/s^2
    C->KAPPA     = C->RADIUS * SQRTST; // m^2/s
    C->DTDP      = 1.0 / (C->DENSITY * C->ENTROPY); // K/Pa
    C->DSDR      = C->ENTROPY / C->RADIUS; // J/kg.K.m
    C->DTDR      = C->TEMP / C->RADIUS; // K/m
    C->GSUPER    = C->GRAVITY * C->DTDR;
    C->VISC      = C->DENSITY * C->KAPPA; // Pa.s
    C->LOG10VISC = PetscLog10Real( C->VISC ); // log10(Pa.s)
    C->COND      = C->ENTROPY * C->DENSITY * C->KAPPA; // W/m.K
    C->SIGMA     = C->FLUX * 1.0 / PetscPowScalar( C->TEMP, 4.0 ); // W/m^2.K^4
    C->LHS       = C->MASS * C->TEMP; // kg.K
    /* TODO: add scaling for internal heat generation since this might
       be required in the internal heat generation functions */
    /* the augmented rhs vector contains various quantities
       with different units, so we cannot scale simply by multiplying
       by a constant value */
    C->RHS       = 1.0; // no scaling, as per comment above
    /* TODO: allow user-specification?  For units of wt % this must
       be 1E2 and for units of ppm this must be 1E6 */
    //C->VOLSCALE = 1.0; // mass fraction
    //C->VOLSCALE  = 1.0E2; // wt %
    //C->VOLSCALE  = 1.0E6; // ppm
    C->VOLSCALE = 1.0E3;

    PetscFunctionReturn(0);
}

/* Initialize Constants, checking for command line arguments.
   For now we only allow a single flag to set all scaling to 1 (hence running
   in "dimensional mode") */
static PetscErrorCode InitializeConstantsAndSetFromOptions(Constants *C)
{
  PetscErrorCode ierr;
  PetscBool      dimensionalMode = PETSC_FALSE;

  PetscFunctionBeginUser;
  ierr = PetscOptionsGetBool(NULL,NULL,"-dimensional",&dimensionalMode,NULL);CHKERRQ(ierr);
  if (dimensionalMode) {
    ierr = SetConstants(C,1.0,1.0,1.0,1.0);CHKERRQ(ierr);
  } else {
    PetscScalar RADIUS0 = 6371000.0; // m
    ierr = PetscOptionsGetScalar(NULL,NULL,"-radius0",&RADIUS0,NULL);CHKERRQ(ierr);  
    PetscScalar ENTROPY0 = 2993.025100070677; // J/kg K
    ierr = PetscOptionsGetScalar(NULL,NULL,"-entropy0",&ENTROPY0,NULL);CHKERRQ(ierr);  
    PetscScalar TEMPERATURE0 = 4033.6070755893948; // K
    ierr = PetscOptionsGetScalar(NULL,NULL,"-temperature0",&TEMPERATURE0,NULL);CHKERRQ(ierr);
    PetscScalar DENSITY0 = 4613.109568155063; // kg/m^3
    ierr = PetscOptionsGetScalar(NULL,NULL,"-density0",&DENSITY0,NULL);CHKERRQ(ierr);
    ierr = SetConstants(C,RADIUS0,TEMPERATURE0,ENTROPY0,DENSITY0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
This function (and subfunctions) should be the only place that
custom command-line parameters (those not defined by PETSc) should be accessed.

Parameters are specified by the user in dimensional (unscaled) form,
but they are all stored in non-dimensional (scaled) form.

 */
PetscErrorCode InitializeParametersAndSetFromOptions(Parameters *P)
{
  PetscErrorCode       ierr;
  AtmosphereParameters *Ap = &P->atmosphere_parameters;
  VolatileParameters   *H2O = &Ap->H2O_volatile_parameters;
  VolatileParameters   *CO2 = &Ap->CO2_volatile_parameters;
  Constants const      *C  = &P->constants;

  PetscFunctionBegin;

  /* Constants (scalings) must be set first, as they are used to scale
     other parameters */
  ierr = InitializeConstantsAndSetFromOptions(&P->constants);CHKERRQ(ierr);

  /* For each entry in parameters, we set a default value and immediately scale it.
     Dimensional/unscaled quantities are not explicitly stored.
     All SI units unless non-dimensional.
     */

  /* Time frame parameters */
  P->maxsteps    = 100000000; /* Effectively infinite */
  P->t0          = 0.0;
  // TODO -dtmacro, contrary to the pattern, are still input in nondim units..
  // TODO: consider *only* accepting -dtmacro_years
  {
    PetscReal dtmacro_years;
    PetscBool dtmacro_set = PETSC_FALSE, dtmacro_years_set = PETSC_FALSE, nstepsmacro_set = PETSC_FALSE,
              early = PETSC_FALSE, middle=PETSC_FALSE, late = PETSC_FALSE;
    P->nstepsmacro = 1000;
    ierr = PetscOptionsGetInt(NULL,NULL,"-nstepsmacro",&P->nstepsmacro,&nstepsmacro_set);CHKERRQ(ierr);
    P->dtmacro = 1000000.0;
    ierr = PetscOptionsGetReal(NULL,NULL,"-dtmacro",&P->dtmacro,&dtmacro_set);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-dtmacro_years",&dtmacro_years,&dtmacro_years_set);CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL,NULL,"-early",&early,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL,NULL,"-middle",&middle,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL,NULL,"-late",&late,NULL);CHKERRQ(ierr);
    if (early || middle || late ){
      if(dtmacro_set || dtmacro_years_set || nstepsmacro_set){
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -early, -middle, or -late provided. Ignoring -nstepsmacro, -nstepmacro_years, and/or -nstepsmacro\n");CHKERRQ(ierr);
      }
      if ((int)early + (int)middle + (int)late > 1) {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"only one of -early, -middle, or -late may be provided");
      }
      /* all these use dtmacro_years, so must turn on flag to activate
         scaling below */
      dtmacro_years_set = PETSC_TRUE;
      if (early) {
        /* early evolution to 10 kyr */
        P->nstepsmacro = 1000;
        dtmacro_years = 10; // years
      } else if (middle) {
        /* middle evolution to 100 Myr */
        P->nstepsmacro = 10000;
        dtmacro_years = 10000; //years
      } else if (late) {
        /* late evolution to 4.55 Byr */
        P->nstepsmacro = 455;
        dtmacro_years = 1000000000000; // years
      }
    }

    if (dtmacro_set && dtmacro_years_set) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: both -dtmacro and -dtmacro_years provided. Using -dtmacro\n");CHKERRQ(ierr);
    /* this is also now called for -early, -middle, and -late */
    } else if (dtmacro_years_set) {
        P->dtmacro = dtmacro_years / C->TIMEYRS;
    }
  }

  /* Grid parameters */
  P->numpts_b = 100;
  //P->numpts_b= 278
  //P->numpts_b= 372
  //P->numpts_b= 656
  //P->numpts_b= 2939
  //P->numpts_b= 5802
  //P->numpts_b= 11532
  //P->numpts_b= 22996
  //P->numpts_b= 45928
  P->numpts_s = P->numpts_b - 1;
  {
    PetscBool set = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&P->numpts_s,&set);CHKERRQ(ierr);
    if (set) P->numpts_b = P->numpts_s + 1;
  }

  /* Monitor */
  P->monitor = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-monitor",&P->monitor,NULL);CHKERRQ(ierr);

  /* Energy terms to include */
  P->CONDUCTION = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-CONDUCTION",&P->CONDUCTION,NULL);CHKERRQ(ierr);
  P->CONVECTION = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-CONVECTION",&P->CONVECTION,NULL);CHKERRQ(ierr);
  P->MIXING = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-MIXING",&P->MIXING,NULL);CHKERRQ(ierr);
  P->SEPARATION = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-SEPARATION",&P->SEPARATION,NULL);CHKERRQ(ierr);
  P->HRADIO = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-HRADIO",&P->HRADIO,NULL);CHKERRQ(ierr);
  P->HTIDAL = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-HTIDAL",&P->HTIDAL,NULL);CHKERRQ(ierr);

  P->mixing_length = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-mixing_length",&P->mixing_length,NULL);CHKERRQ(ierr);

  /* initial entropy at top of adiabat (J/kg-K) */
  P->sinit = 3052.885602072091;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-sinit",&P->sinit,NULL);CHKERRQ(ierr);
  P->sinit /= C->ENTROPY;

  /* initial entropy gradient (J/kg-K-m) */
  P->ic_dsdr = -4.6978890285209187e-07;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-ic_dsdr",&P->ic_dsdr,NULL);CHKERRQ(ierr);
  P->ic_dsdr /= C->DSDR;

  /* radius of planet (m) */
  P->radius = 6371000.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-radius",&P->radius,NULL);CHKERRQ(ierr);
  P->radius /= C->RADIUS;

  /* core size (non-dimensional) */
  P->coresize = 0.55;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-coresize",&P->coresize,NULL);CHKERRQ(ierr);

  /* surface density (kg/m^3) for Adams-Williamson EOS for pressure */
  P->rhos = 4078.95095544;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-rhos",&P->rhos,NULL);CHKERRQ(ierr);
  P->rhos /= C->DENSITY;

  /* parameter (1/m) for Adams-Williamson EOS for pressure */
  P->beta = 1.1115348931000002e-07;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-beta",&P->beta,NULL);CHKERRQ(ierr);
  P->beta *= C->RADIUS;

  /* grain size (m) */
  P->grain = 1.0E-3;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-grain",&P->grain,NULL);CHKERRQ(ierr);
  P->grain /= C->RADIUS;

  /* gravity (m/s^2), must be negative */
  P->gravity = -10.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-gravity",&P->gravity,NULL);CHKERRQ(ierr);
  P->gravity /= C->GRAVITY;

  /* melt fraction threshold for rheology */
  P->phi_critical = 0.4; // non dimensional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-phi_critical",&P->phi_critical,NULL);CHKERRQ(ierr);
  /* melt fraction transition width for rheology */
  /* when PHI_WIDTH = 0.15, oscillations appear in the volatile contents
     due to a rigid crust forming at the top of the model.  Reducing the value
     to 0.2 helps to alleviate this problem.  So evidently the viscosity contrast
     across nodes matters. */
  P->phi_width = 0.15; // non dimensional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-phi_width",&P->phi_width,NULL);CHKERRQ(ierr);

  /* melt fraction shape transition for skew */
  // FIXME: I think this is broken at present
  P->phi_skew = 0.0; // non dimensional

  /* core density (kg/m^3) */
  P->rho_core = 10738.332568062382;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-rho_core",&P->rho_core,NULL);CHKERRQ(ierr);
  P->rho_core /= C->DENSITY;

  /* heat capacity of core (J/kg-K) */
  P->cp_core = 880.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-cp_core",&P->cp_core,NULL);CHKERRQ(ierr);
  P->cp_core /= C->ENTROPY;

  /* mass-weighted average core temperature as a fraction */
  /* of CMB temperature (non-dimensional) */
  P->tfac_core_avg = 1.147;

  /* smoothing width */
  P->swidth = 1.0E-2;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-swidth",&P->swidth,NULL);CHKERRQ(ierr);

  /* solid viscosity (Pa.s) */
  P->log10visc_sol = 21.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-log10visc_sol",&P->log10visc_sol,NULL);CHKERRQ(ierr);
  P->log10visc_sol -= C->LOG10VISC;

  /* solid conductivity (W/m-K) */
  P->cond_sol = 4.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-cond_sol",&P->cond_sol,NULL);CHKERRQ(ierr);
  P->cond_sol /= C->COND;

  /* melt viscosity (Pa.s) */
  P->log10visc_mel = 2.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-log10visc_mel",&P->log10visc_mel,NULL);CHKERRQ(ierr);
  P->log10visc_mel -= C->LOG10VISC;

  /* melt conductivity (W/m-K) */
  P->cond_mel = 4.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-cond_mel",&P->cond_mel,NULL);CHKERRQ(ierr);
  P->cond_mel /= C->COND;

  /* atmosphere parameters
     MODEL = MO_ATMOSPHERE_TYPE_GREY_BODY: grey-body
     MODEL = MO_ATMOSPHERE_TYPE_ZAHNLE: steam atmosphere
     MODEL = MO_ATMOSPHERE_TYPE_VOLATILES: self-consistent atmosphere evolution
     with CO2 and H2O volatiles
     uses plane-parallel radiative eqm model
     of Abe and Matsui (1985)
     */
  Ap->MODEL=MO_ATMOSPHERE_TYPE_GREY_BODY;
  {
    PetscInt  MODEL = 0;
    PetscBool MODELset = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL,NULL,"-MODEL",&MODEL,&MODELset);CHKERRQ(ierr);
    if( MODELset ) Ap->MODEL = MODEL;
  }

  Ap->HYBRID = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-HYBRID",&Ap->HYBRID,NULL);CHKERRQ(ierr);

  /* emissivity is constant for MODEL != MO_ATMOSPHERE_TYPE_VOLATILES */
  Ap->emissivity0 = 1.0; // non-dimensional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-emissivity0",&Ap->emissivity0,NULL);CHKERRQ(ierr);

  /* Stefan-Boltzmann constant (W/m^2K^4) */
  Ap->sigma = 5.670367e-08;
  Ap->sigma /= C->SIGMA;

  /* equilibrium temperature of the planet (K) */
  Ap->teqm = 273.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-teqm",&Ap->teqm,NULL);CHKERRQ(ierr);
  Ap->teqm /= C->TEMP;

  /* for radiative boundary condition at the top surface
     dT = param_utbl_const * [Surface temperature]**3 */
  Ap->PARAM_UTBL = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-PARAM_UTBL",&Ap->PARAM_UTBL,NULL);CHKERRQ(ierr);
  if (Ap->PARAM_UTBL){
      Ap->param_utbl_const = 1.0E-7;
      ierr = PetscOptionsGetScalar(NULL,NULL,"-param_utbl_const",&Ap->param_utbl_const,NULL);CHKERRQ(ierr);
      Ap->param_utbl_const *= PetscSqr(C->TEMP);
  }
  else{
      Ap->param_utbl_const = 0.0;
  }

  Ap->SOLVE_FOR_VOLATILES = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-SOLVE_FOR_VOLATILES",&Ap->SOLVE_FOR_VOLATILES,NULL);CHKERRQ(ierr);

  /* below here are only used for MODEL = MO_ATMOSPHERE_TYPE_VOLATILES */

  /* atmosphere reference pressure (Pa) */
  Ap->P0 = 101325.0; // Pa (= 1 atm)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-P0",&Ap->P0,NULL);CHKERRQ(ierr);
  Ap->P0 /= C->PRESSURE;

  /* H2O volatile */
  H2O->initial = 0.0; // units according to VOLSCALE (typically wt % or ppm)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-H2O_initial",&H2O->initial,NULL);CHKERRQ(ierr);
  H2O->kdist = 1.0E-4; // non-dimensional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-H2O_kdist",&H2O->kdist,NULL);CHKERRQ(ierr);
  // TODO: water saturation limit of 10 ppm?
  H2O->kabs = 0.01; // m^2/kg
  ierr = PetscOptionsGetScalar(NULL,NULL,"-H2O_kabs",&H2O->kabs,NULL);CHKERRQ(ierr);
  H2O->kabs *= C->DENSITY * C->RADIUS;
  H2O->henry = 6.8E-8; // must be mass fraction/Pa
  ierr = PetscOptionsGetScalar(NULL,NULL,"-H2O_henry",&H2O->henry,NULL);CHKERRQ(ierr);
  /* scaled henry constant used in code */
  H2O->henry *= C->VOLSCALE;
  H2O->henry_pow = 1.4285714285714286; // (1.0/0.7)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-H2O_henry_pow",&H2O->henry_pow,NULL);CHKERRQ(ierr);

  /* CO2 volatile */
  CO2->initial = 0.0; // units according to VOLSCALE
  ierr = PetscOptionsGetScalar(NULL,NULL,"-CO2_initial",&CO2->initial,NULL);CHKERRQ(ierr);
  CO2->kdist = 5.0E-4; // non-dimensional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-CO2_kdist",&CO2->kdist,NULL);CHKERRQ(ierr);
  // TODO: water saturation limit of 0.03 ppm
  CO2->kabs = 0.05; // m^2/kg
  ierr = PetscOptionsGetScalar(NULL,NULL,"-CO2_kabs",&CO2->kabs,NULL);CHKERRQ(ierr);
  CO2->kabs *= C->DENSITY * C->RADIUS;
  CO2->henry = 4.4E-12; // must be mass fraction/Pa
  ierr = PetscOptionsGetScalar(NULL,NULL,"-CO2_henry",&CO2->henry,NULL);CHKERRQ(ierr);
  /* scaled henry constant used in code */
  CO2->henry *= C->VOLSCALE;
  CO2->henry_pow = 1.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-CO2_henry_pow",&CO2->henry_pow,NULL);CHKERRQ(ierr);

  /* Get lookup tables */
  ierr = SetLookups(P);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PrintParameters(Parameters const *P)
{
  PetscErrorCode ierr;
  Constants const *C = &P->constants;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"**************** Magma Ocean | Parameters **************\n\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15s %s\n"                 ,"[Scaling]"  ,"","Value"  ,"Units"            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Radius"     ,"",C->RADIUS,"m"                );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Temperature","",C->TEMP  ,"K"                );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Entropy"    ,"",C->ENTROPY,"J/kg-K"          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Density"    ,"",C->DENSITY,"kg/m^3"          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s (%.6g years)\n","Time"       ,"",C->TIME  ,"s",C->TIMEYRS     );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Area"       ,"",C->AREA  ,"m^2"              );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Volume"     ,"",C->VOLUME,"m^3");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Mass"       ,"",C->MASS,"kg");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"SEnergy"    ,"",C->SENERGY,"J/kg");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Energy"     ,"",C->ENERGY,"J");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Pressure"   ,"",C->PRESSURE,"Pa");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Power"      ,"",C->POWER,"W");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Flux"       ,"",C->FLUX,"W/m^2");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dP/dR"      ,"",C->DPDR,"Pa/m");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Gravity"    ,"",C->GRAVITY,"m/s^2");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Kappa"      ,"",C->KAPPA,"m^2/s");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dT/dP"      ,"",C->DTDP,"K/Pa");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dS/dR"      ,"",C->DSDR,"J/kg-K-m");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dT/dR"      ,"",C->DTDR,"K/m");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"GSuper"     ,"",C->GSUPER,"K/s^2");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Visc"       ,"",C->VISC,"Pa-s");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Log10 Visc" ,"",C->LOG10VISC,"log10(Pa-s)");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Cond"       ,"",C->COND,"W/m-K");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Sigma"      ,"",C->SIGMA,"W/m^2-K^4");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Lhs"        ,"",C->LHS,"kg-K");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Rhs"        ,"",C->RHS,"W/kg-K");CHKERRQ(ierr);
  // TODO finish adding units, do other cleanup

  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"                                                                                    );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15s %s\n"    ,"[Parameter]","Non-dim. Value","Dim. Value"       ,"Units" );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15.6g %-15.6g %s\n","dtmacro"    ,P->dtmacro      ,P->dtmacro*C->TIME ,"s"     );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15d\n"             ,"nstepsmacro",P->nstepsmacro                               );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15d\n"             ,"numpts_b"   ,P->numpts_b                                  );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15d\n"             ,"numpts_s"   ,P->numpts_s                                  );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15.6g %-15.6g %s\n","S_init"     ,P->sinit        ,P->sinit*C->ENTROPY,"J/kg-K");CHKERRQ(ierr);
  // TODO add the rest..
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"liquidus data file"          ,P->liquidusFilename         );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"solidus data file"           ,P->solidusFilename          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"alphaSol data file"          ,P->alphaSolFilename         );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"alphaMel data file"          ,P->alphaMelFilename         );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"cpSol data file"             ,P->cpSolFilename            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"cpMel data file"             ,P->cpMelFilename            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"dtdpsSol data file"          ,P->dtdpsSolFilename         );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"dtdpsMel data file"          ,P->dtdpsMelFilename         );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"rhoSol data file"            ,P->rhoSolFilename           );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"rhoMel data file"            ,P->rhoMelFilename           );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"tempSol data file"           ,P->tempSolFilename          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"tempMel data file"           ,P->tempMelFilename          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************************\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Helper routine to prepend the root directory to a relative path */
static PetscErrorCode MakeRelativePathAbsolute(char* path) {
  PetscErrorCode ierr;
  char tmp[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  ierr = PetscStrcpy(tmp,path);CHKERRQ(ierr);
  ierr = PetscStrcpy(path,MAGMA_ROOT_DIR_STR);CHKERRQ(ierr); /* lookup.h */
  ierr = PetscStrcat(path,"/");CHKERRQ(ierr); /* not portable */
  ierr = PetscStrcat(path,tmp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SetLookups( Parameters *P )
{
    PetscErrorCode ierr;
    Constants *C = &P->constants;

    PetscFunctionBeginUser;
    /* set all 1-D and 2-D lookups */
#if (defined VERBOSE)
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"set_lookup:\n" );CHKERRQ(ierr);
    }
#endif

    /* Based on a user-supplied flag, determine which files to load.
       Note that we prepend a string, MAGMA_ROOT_DIR_STR, and /, to make the
       paths absolute, so the code can be run (tested) from anywhere.
       See lookup.h for the paths */
    {
      char lscTypeString[PETSC_MAX_PATH_LEN] = "default";
      PetscBool isDefault=PETSC_FALSE,isStixrude2009=PETSC_FALSE,isAndrault2011=PETSC_FALSE;
      ierr = PetscOptionsGetString(NULL,NULL,"-curves",lscTypeString,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
      ierr = PetscStrcmp(lscTypeString,"default",&isDefault);CHKERRQ(ierr);
      ierr = PetscStrcmp(lscTypeString,"stixrude2009",&isStixrude2009);CHKERRQ(ierr);
      ierr = PetscStrcmp(lscTypeString,"andrault2011",&isAndrault2011);CHKERRQ(ierr);
      if (!(isDefault || isStixrude2009 || isAndrault2011)){
        ierr = PetscPrintf(PETSC_COMM_WORLD,"*************** WARNING ***************\n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Unrecognized -curves choice %s provided. Using defaults.\nCurrent options include:\n -curves stixrude2009\n -curves andrault2011\n",lscTypeString);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"***************************************\n");CHKERRQ(ierr);
      }
      if (isStixrude2009) {
        ierr = PetscStrcpy(P->liquidusFilename,LIQUIDUS_STIXRUDE2009);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(P->liquidusFilename);CHKERRQ(ierr);
        ierr = PetscStrcpy(P->solidusFilename,SOLIDUS_STIXRUDE2009);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(P->solidusFilename);CHKERRQ(ierr);
      } else if (isAndrault2011 || isDefault) { /* Default */
        ierr = PetscStrcpy(P->liquidusFilename,LIQUIDUS_ANDRAULT2011);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(P->liquidusFilename);CHKERRQ(ierr);
        ierr = PetscStrcpy(P->solidusFilename,SOLIDUS_ANDRAULT2011);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(P->solidusFilename);CHKERRQ(ierr);
      } else {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-curves logic processing failed");
      }

      ierr = PetscStrcpy(P->alphaSolFilename,ALPHA_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->alphaSolFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(P->cpSolFilename,CP_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->cpSolFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(P->dtdpsSolFilename,DTDPS_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->dtdpsSolFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(P->rhoSolFilename,RHO_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->rhoSolFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(P->tempSolFilename,TEMP_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->tempSolFilename);CHKERRQ(ierr);

      ierr = PetscStrcpy(P->alphaMelFilename,ALPHA_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->alphaMelFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(P->cpMelFilename,CP_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->cpMelFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(P->dtdpsMelFilename,DTDPS_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->dtdpsMelFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(P->rhoMelFilename,RHO_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->rhoMelFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(P->tempMelFilename,TEMP_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(P->tempMelFilename);CHKERRQ(ierr);
    }

    /* solid lookups */
    /* 2d */
    set_interp2d( P->alphaSolFilename, &P->solid_prop.alpha, C->PRESSURE, C->ENTROPY, 1.0/C->TEMP );
    set_interp2d( P->cpSolFilename, &P->solid_prop.cp, C->PRESSURE, C->ENTROPY, C->ENTROPY );
    set_interp2d( P->dtdpsSolFilename, &P->solid_prop.dTdPs, C->PRESSURE, C->ENTROPY, C->DTDP );
    set_interp2d( P->rhoSolFilename, &P->solid_prop.rho, C->PRESSURE, C->ENTROPY, C->DENSITY );
    set_interp2d( P->tempSolFilename, &P->solid_prop.temp, C->PRESSURE, C->ENTROPY, C->TEMP );
    /* const */
    // FIXME: remove this redundancy
    P->solid_prop.cond = P->cond_sol;
    P->solid_prop.log10visc = P->log10visc_sol;

    /* melt lookups */
    /* 2d */
    set_interp2d( P->alphaMelFilename, &P->melt_prop.alpha, C->PRESSURE, C->ENTROPY, 1.0/C->TEMP );
    set_interp2d( P->cpMelFilename, &P->melt_prop.cp, C->PRESSURE, C->ENTROPY, C->ENTROPY );
    set_interp2d( P->dtdpsMelFilename, &P->melt_prop.dTdPs, C->PRESSURE, C->ENTROPY, C->DTDP );
    set_interp2d( P->rhoMelFilename, &P->melt_prop.rho, C->PRESSURE, C->ENTROPY, C->DENSITY );
    set_interp2d( P->tempMelFilename, &P->melt_prop.temp, C->PRESSURE, C->ENTROPY, C->TEMP );
    /* const */
    // FIXME remove this redundancy
    P->melt_prop.cond = P->cond_mel;
    P->melt_prop.log10visc = P->log10visc_mel;

    /* liquidus and solidus */
    /* 1d */

    set_interp1d( P->liquidusFilename, &P->solid_prop.liquidus, NLS, C->PRESSURE, C->ENTROPY );
    set_interp1d( P->liquidusFilename, &P->melt_prop.liquidus, NLS, C->PRESSURE, C->ENTROPY );
    /* duplication here, but want to remain flexible for future
       approaches for a multicomponent system */
    set_interp1d( P->solidusFilename, &P->solid_prop.solidus, NLS, C->PRESSURE, C->ENTROPY );
    set_interp1d( P->solidusFilename, &P->melt_prop.solidus, NLS, C->PRESSURE, C->ENTROPY );

    PetscFunctionReturn(0);
}
