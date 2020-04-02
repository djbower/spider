/*
Parameter Management

Parameters should only ever be set by the functions in this file. That is, everywhere else they should be considered read-only.

Custom PETSc command line options should only ever be parsed here.
*/

#include "parameters.h"
#include "ctx.h"
#include "ic.h"
#include "eos.h"
// FIXME
//#include "composition.h"

static PetscErrorCode set_start_time_from_file( Parameters * , const char * );


static PetscErrorCode SetConstants( Constants *C, PetscReal RADIUS, PetscReal TEMPERATURE, PetscReal ENTROPY, PetscReal DENSITY, PetscReal VOLATILE )
{
    PetscScalar SQRTST;

    PetscFunctionBeginUser;
    /* 29 constants to set (excluding SQRTST which is a convenience
       parameter) */
    SQRTST = PetscSqrtScalar( ENTROPY * TEMPERATURE );

    C->RADIUS    = RADIUS; // m
    C->TEMP      = TEMPERATURE; // K
    C->ENTROPY   = ENTROPY; // (specific) J/kg.K
    C->DENSITY   = DENSITY; // kg/m^3
    C->VOLATILE  = VOLATILE;
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
    C->GRAVITY   = C->ENTROPY * C->TEMP / C->RADIUS; // m/s^2
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
    C->HEATGEN   = PetscPowScalar( C->ENTROPY*C->TEMP, 3.0/2.0 ) / C->RADIUS; // W/kg
    /* the full rhs vector contains various quantities
       with different units, so we cannot scale simply by multiplying
       by a constant value */
    C->RHS       = 1.0; // no scaling, as per comment above

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
    ierr = SetConstants(C,1.0,1.0,1.0,1.0,1.0);CHKERRQ(ierr);
  } else {
    PetscScalar RADIUS0 = 6371000.0; // m
    ierr = PetscOptionsGetScalar(NULL,NULL,"-radius0",&RADIUS0,NULL);CHKERRQ(ierr);
    PetscScalar ENTROPY0 = 2993.025100070677; // J/kg K
    ierr = PetscOptionsGetScalar(NULL,NULL,"-entropy0",&ENTROPY0,NULL);CHKERRQ(ierr);
    PetscScalar TEMPERATURE0 = 4033.6070755893948; // K
    ierr = PetscOptionsGetScalar(NULL,NULL,"-temperature0",&TEMPERATURE0,NULL);CHKERRQ(ierr);
    PetscScalar DENSITY0 = 4613.109568155063; // kg/m^3
    ierr = PetscOptionsGetScalar(NULL,NULL,"-density0",&DENSITY0,NULL);CHKERRQ(ierr);
    PetscScalar VOLATILE0 = 1.0;
    ierr = PetscOptionsGetScalar(NULL,NULL,"-volatile0",&VOLATILE0,NULL);CHKERRQ(ierr);
    ierr = SetConstants(C,RADIUS0,TEMPERATURE0,ENTROPY0,DENSITY0,VOLATILE0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode VolatileParametersSetFromOptions(VolatileParameters *vp, const Constants *C)
{
  PetscErrorCode ierr;
  char           buf[1024]; /* max size */
  PetscBool      set;

  PetscFunctionBeginUser;
  /* Accept -prefix_YYY to populate vp->YYY. Most are required and an error is thrown
     if they are missing. Note that this code has a lot of duplication */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_initial_total_abundance");CHKERRQ(ierr);
  vp->initial_total_abundance = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf, &vp->initial_total_abundance,&set);CHKERRQ(ierr);
  vp->initial_total_abundance /= 1.0E6 * C->VOLATILE;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_initial_atmos_pressure");CHKERRQ(ierr);
  vp->initial_atmos_pressure = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf, &vp->initial_atmos_pressure,&set);CHKERRQ(ierr);
  vp->initial_atmos_pressure /= C->PRESSURE;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_kdist");CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->kdist,&set);CHKERRQ(ierr);
  if (!set) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing argument %s",buf);

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_kabs");CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->kabs,&set);CHKERRQ(ierr);
  if (!set) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing argument %s",buf);
  vp->kabs *= C->DENSITY * C->RADIUS;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_SOLUBILITY");CHKERRQ(ierr);
  vp->SOLUBILITY = 1; // Modified Henry's law
  ierr = PetscOptionsGetInt(NULL,NULL,buf,&vp->SOLUBILITY,&set);CHKERRQ(ierr);

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_henry_pow");CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->henry_pow,&set);CHKERRQ(ierr);
  if (!set) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing argument %s",buf);
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_henry");CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->henry,&set);CHKERRQ(ierr);
  if (!set) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing argument %s",buf);
  /* ensure that V->dxdp calculation does not return NaN when vp->henry is zero */
  if(vp->henry == 0){
    vp->henry_pow = 1.0;
  }
  vp->henry /= 1.0E6 * C->VOLATILE * PetscPowScalar(C->PRESSURE, -1.0/vp->henry_pow);

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_jeans_value");CHKERRQ(ierr);
  vp->jeans_value = 0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->jeans_value,&set);CHKERRQ(ierr);

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_R_thermal_escape_value");CHKERRQ(ierr);
  vp->R_thermal_escape_value = 0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->R_thermal_escape_value,&set);CHKERRQ(ierr);

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_constant_escape_value");CHKERRQ(ierr);
  vp->constant_escape_value = 0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->constant_escape_value,&set);CHKERRQ(ierr);
  vp->constant_escape_value *= C->TIME / (C->VOLATILE * C->MASS);
  vp->constant_escape_value /= 4.0 * PETSC_PI;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_molar_mass");CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->molar_mass,&set);CHKERRQ(ierr);
  if (!set) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing argument %s",buf);
  vp->molar_mass /= C->MASS;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_cross_section");CHKERRQ(ierr);
  vp->cross_section = 1.0E-18; // m^2, Johnson et al. (2015), N2+N2 collisions
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->cross_section,NULL);CHKERRQ(ierr);
  vp->cross_section /= C->AREA;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",vp->prefix,"_poststep_change");CHKERRQ(ierr);
  vp->poststep_change = -1; // fractional (negative value is OFF)
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&vp->poststep_change,&set);CHKERRQ(ierr);
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
  Constants const      *C  = &P->constants;
  RadiogenicIsotopeParameters *al26 = &P->al26_parameters;
  RadiogenicIsotopeParameters *k40 = &P->k40_parameters;
  RadiogenicIsotopeParameters *fe60 = &P->fe60_parameters;
  RadiogenicIsotopeParameters *th232 = &P->th232_parameters;
  RadiogenicIsotopeParameters *u235 = &P->u235_parameters;
  RadiogenicIsotopeParameters *u238 = &P->u238_parameters;
  PetscInt             v;

  // FIXME
  //CompositionParameters      *Compp = &P->composition_parameters;

  PetscFunctionBegin;

  /* Constants (scalings) must be set first, as they are used to scale
     other parameters */
  ierr = InitializeConstantsAndSetFromOptions(&P->constants);CHKERRQ(ierr);

  /* Must set EOS after setting constants, but before boundary conditions
     since EOS might be required to map temperature to entropy
     and vice versa for boundary conditions */
  /* TODO: if we have an option of analytical EOSs, we might not need to set the lookups
     for every model run */

  P->SOLID_EOS = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-SOLID_EOS",&P->SOLID_EOS,NULL);CHKERRQ(ierr);

  P->MELT_EOS = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-MELT_EOS",&P->MELT_EOS,NULL);CHKERRQ(ierr);

  /* function is in eos.c */
  ierr = set_eos( P );CHKERRQ(ierr);

  /* For each entry in parameters, we set a default value and immediately scale it.
     Dimensional/unscaled quantities are not explicitly stored.
     All SI units unless non-dimensional.
     */

  /* Time frame parameters */
  P->maxsteps    = 100000000; /* Effectively infinite */

  P->nstepsmacro = 18;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nstepsmacro",&P->nstepsmacro,NULL);CHKERRQ(ierr);

  /* start time (years) P->t0 is set further down, since it may
     acquire a value from a restart value */

  /* step time (years) */
  P->dtmacro = 100;
  ierr = PetscOptionsGetReal(NULL,NULL,"-dtmacro",&P->dtmacro,NULL);CHKERRQ(ierr);
  P->dtmacro /= C->TIMEYRS; // non-dimensional for time stepping

  /* Grid parameters */
  P->numpts_b = 200;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&P->numpts_b,NULL);CHKERRQ(ierr);
  P->numpts_s = P->numpts_b - 1;

  /* RollBack and PostStep options */
  P->rollBackActive = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-activate_rollback",&P->rollBackActive,NULL);CHKERRQ(ierr);
  P->postStepActive = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-activate_poststep",&P->postStepActive,NULL);CHKERRQ(ierr);

  /* Output Options */
  P->monitor = PETSC_TRUE;
  ierr = PetscStrcpy(P->outputDirectory,"output");CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-outputDirectory",P->outputDirectory,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);

  P->SOLID_CONVECTION_ONLY = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-SOLID_CONVECTION_ONLY",&P->SOLID_CONVECTION_ONLY,NULL);CHKERRQ(ierr);

  P->LIQUID_CONVECTION_ONLY = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-LIQUID_CONVECTION_ONLY",&P->LIQUID_CONVECTION_ONLY,NULL);CHKERRQ(ierr);

  /* Energy terms to include */
  P->CONDUCTION = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-CONDUCTION",&P->CONDUCTION,NULL);CHKERRQ(ierr);
  P->CONVECTION = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-CONVECTION",&P->CONVECTION,NULL);CHKERRQ(ierr);
  P->HRADIO = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-HRADIO",&P->HRADIO,NULL);CHKERRQ(ierr);
  P->HTIDAL = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-HTIDAL",&P->HTIDAL,NULL);CHKERRQ(ierr);

  P->MIXING = PETSC_TRUE;
  P->SEPARATION = PETSC_TRUE;
  /* mixed phase energy terms */
  if( P->SOLID_CONVECTION_ONLY || P->LIQUID_CONVECTION_ONLY ){
      P->MIXING = PETSC_FALSE;
      P->SEPARATION = PETSC_FALSE;
  }
  ierr = PetscOptionsGetBool(NULL,NULL,"-MIXING",&P->MIXING,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-SEPARATION",&P->SEPARATION,NULL);CHKERRQ(ierr);

  P->mixing_length = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-mixing_length",&P->mixing_length,NULL);CHKERRQ(ierr);
  P->mixing_length_layer_radius = 0.0;
  if ( P->mixing_length==3 ){
    ierr = PetscOptionsGetScalar(NULL,NULL,"-mixing_length_layer_radius",&P->mixing_length_layer_radius,NULL);CHKERRQ(ierr);
  }

  P->Mg_Si0 = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-Mg_Si0",&P->Mg_Si0,NULL);CHKERRQ(ierr);
  P->Mg_Si1 = 0.0;
  if ( P->mixing_length == 3){
    ierr = PetscOptionsGetScalar(NULL,NULL,"-Mg_Si1",&P->Mg_Si1,NULL);CHKERRQ(ierr);
  }

  /* start time (years) */
  P->t0 = 0.0;

  /* initial condition for interior */
  P->IC_INTERIOR = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-IC_INTERIOR",&P->IC_INTERIOR,NULL);CHKERRQ(ierr);

  ierr = PetscStrcpy(P->ic_interior_filename,"restart.json"); CHKERRQ(ierr);
  if ( P->IC_INTERIOR==2 ){
    ierr = PetscOptionsGetString(NULL,NULL,"-ic_interior_filename",P->ic_interior_filename,PETSC_MAX_PATH_LEN,NULL); CHKERRQ(ierr);
    /* set start time from restart file */
    /* TODO: get time from interior restart file, but could add
       other options to get time from e.g. atmosphere restart file */
    ierr = set_start_time_from_file( P, P->ic_interior_filename ); CHKERRQ(ierr);
  }
  else{
    PetscScalar   t0 = 0.0;
    PetscBool t0_set = PETSC_FALSE;
    ierr = PetscOptionsGetReal(NULL,NULL,"-t0",&t0,&t0_set);CHKERRQ(ierr);
    if( t0_set ) P->t0 = t0;
    P->t0 /= C->TIMEYRS; // non-dimensional for time stepping
  }

  P->ic_melt_pressure = 30.0; // GPa
  if ( P->IC_INTERIOR==3 ){
    ierr = PetscOptionsGetScalar(NULL,NULL,"-ic_melt_pressure",&P->ic_melt_pressure,NULL); CHKERRQ(ierr);
  }

  /* initial entropy at top of adiabat (J/kg-K) */
  P->ic_adiabat_entropy = 3052.885602072091;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-ic_adiabat_entropy",&P->ic_adiabat_entropy,NULL);CHKERRQ(ierr);
  P->ic_adiabat_entropy /= C->ENTROPY;

  /* initial entropy gradient (J/kg-K-m) */
  P->ic_dsdr = -4.6978890285209187e-07;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-ic_dsdr",&P->ic_dsdr,NULL);CHKERRQ(ierr);
  P->ic_dsdr /= C->DSDR;

  /* initial entropy at the surface, or set to P->ic_adiabat_entropy if P->ic_surface_entropy < 0.0 */
  P->ic_surface_entropy = -1;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-ic_surface_entropy",&P->ic_surface_entropy,NULL);CHKERRQ(ierr);
  if( P->ic_surface_entropy > 0.0 ){
    P->ic_surface_entropy /= C->ENTROPY;
  }

  /* initial entropy at the core-mantle boundary, or leave as value at base of adiabat if P->ic_core_entropy < 0.0 */
  P->ic_core_entropy = -1;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-ic_core_entropy",&P->ic_core_entropy,NULL);CHKERRQ(ierr);
  if( P->ic_core_entropy > 0.0 ){
    P->ic_core_entropy /= C->ENTROPY;
  }

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

  /* core mass (calculated, non-dimensional scaled mass) */
  P->coremass = 1.0/3.0 * PetscPowScalar( P->coresize, 3.0) * PetscPowScalar( P->radius, 3.0 );
  P->coremass *= P->rho_core;

  /* heat capacity of core (J/kg-K) */
  P->cp_core = 880.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-cp_core",&P->cp_core,NULL);CHKERRQ(ierr);
  P->cp_core /= C->ENTROPY;

  /* mass-weighted average core temperature as a fraction */
  /* of CMB temperature (non-dimensional) */
  P->tfac_core_avg = 1.147;

  /* smoothing width */
  P->matprop_smooth_width = 1.0E-2;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-matprop_smooth_width",&P->matprop_smooth_width,NULL);CHKERRQ(ierr);

  /* TODO: primary viscosity is stored with the EOS parameters, and may these should be too.  Or setup
     a viscosity struct to contain viscosity-related parameters */

  /* solid activation energy (J/mol) */
  PetscScalar Rgas = 8.314; // gas constant (J/K/mol)
  P->activation_energy_sol = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-activation_energy_sol",&P->activation_energy_sol,NULL);CHKERRQ(ierr);
  P->activation_energy_sol /= Rgas * C->TEMP; // this is a new energy scale (i.e., not C->ENERGY defined above)

  /* solid activation volume (m^3/mol) */
  /* The numerical value in units of m^3/mol is the same as that in units of J/mol/Pa */
  /* You can convince yourself of this by using the scalings for ENERGY and PRESSURE to
     see that this is true */
  P->activation_volume_sol = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-activation_volume_sol",&P->activation_volume_sol,NULL);CHKERRQ(ierr);
  P->activation_volume_sol *= C->PRESSURE;
  P->activation_volume_sol /= Rgas * C->TEMP;

  /* viscosity cut-offs */
  P->log10visc_min = -1.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-log10visc_min",&P->log10visc_min,NULL);CHKERRQ(ierr);
  if( P->log10visc_min > 0.0 ){
      P->log10visc_min -= C->LOG10VISC;
  }
  P->log10visc_max = -1.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-log10visc_max",&P->log10visc_max,NULL);CHKERRQ(ierr);
  if( P->log10visc_max > 0.0 ){
      P->log10visc_max -= C->LOG10VISC;
  }

  /* option to scale eddy diffusivities for temperature and chemistry, or set as constants */
  P->eddy_diffusivity_thermal = 1.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-eddy_diffusivity_thermal",&P->eddy_diffusivity_thermal,NULL);CHKERRQ(ierr);
  /* if input is negative, then set as constant.  Retain negative sign to use as flag */
  if( P->eddy_diffusivity_thermal < 0.0){
    /* must scale */
    P->eddy_diffusivity_thermal /= C->KAPPA;
  }
  /* otherwise, we just scale the calculated eddy diffusivity by the user-specified constant */

  P->eddy_diffusivity_chemical = 1.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-eddy_diffusivity_chemical",&P->eddy_diffusivity_chemical,NULL);CHKERRQ(ierr);
  /* if input is negative, then set as constant.  Retain negative sign to use as flag */
  if( P->eddy_diffusivity_chemical < 0.0){
    /* must scale */
    P->eddy_diffusivity_chemical /= C->KAPPA;
  }
  /* otherwise, we just scale the calculated eddy diffusivity by the user-specified constant */

  /* viscous lid added by Rob Spaargaren */
  P->VISCOUS_LID = 0;
  ierr = PetscOptionsGetInt(NULL,NULL,"-VISCOUS_LID",&P->VISCOUS_LID,NULL);CHKERRQ(ierr);
  P->lid_log10visc = 0.0;
  P->lid_thickness = 0.0; /* metres */
  if ( P->VISCOUS_LID ){
      ierr = PetscOptionsGetScalar(NULL,NULL,"-lid_log10visc",&P->lid_log10visc,NULL);CHKERRQ(ierr);
      ierr = PetscOptionsGetScalar(NULL,NULL,"-lid_thickness",&P->lid_thickness,NULL);CHKERRQ(ierr);
      P->lid_thickness /= C->RADIUS;
  }

  /* core boundary condition */
  P->CORE_BC=MO_CORE_TYPE_COOLING;
  {
    PetscInt  CORE_BC = 0;
    PetscBool CORE_BCset = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL,NULL,"-CORE_BC",&CORE_BC,&CORE_BCset);CHKERRQ(ierr);
    if( CORE_BCset ) P->CORE_BC = CORE_BC;
  }
  P->core_bc_value = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-core_bc_value",&P->core_bc_value,NULL);CHKERRQ(ierr);
  switch( P->CORE_BC ){
    case 1:
      // CORE_BC = MO_CORE_TYPE_COOLING: simple core cooling
      P->core_bc_value = 0.0;
      break;
    case 2:
      // CORE_BC = MO_CORE_TYPE_HEAT_FLUX: heat flux (prescribed)
      P->core_bc_value /= C->FLUX;
      break;
    case 3:
      /* CORE_BC = MO_CORE_TYPE_ENTROPY: entropy
         do nothing */
      break;
    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported CORE_BC value %d provided",P->CORE_BC);
      break;
  }

  /* initial condition for atmosphere */
  Ap->IC_ATMOSPHERE = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-IC_ATMOSPHERE",&Ap->IC_ATMOSPHERE,NULL);CHKERRQ(ierr);

  ierr = PetscStrcpy(Ap->ic_atmosphere_filename,"restart.json"); CHKERRQ(ierr);
  if ( (Ap->IC_ATMOSPHERE==2) || (Ap->IC_ATMOSPHERE==3) ){
    ierr = PetscOptionsGetString(NULL,NULL,"-ic_atmosphere_filename",Ap->ic_atmosphere_filename,PETSC_MAX_PATH_LEN,NULL); CHKERRQ(ierr);
  }

  /* (top) surface boundary condition */
  Ap->SURFACE_BC=MO_ATMOSPHERE_TYPE_GREY_BODY;
  {
    PetscInt  SURFACE_BC = 0;
    PetscBool SURFACE_BCset = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL,NULL,"-SURFACE_BC",&SURFACE_BC,&SURFACE_BCset);CHKERRQ(ierr);
    if( SURFACE_BCset ) Ap->SURFACE_BC = SURFACE_BC;
  }
  Ap->surface_bc_value = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-surface_bc_value",&Ap->surface_bc_value,NULL);CHKERRQ(ierr);
  switch( Ap->SURFACE_BC ){
    case 1:
      /* SURFACE_BC = MO_ATMOSPHERE_TYPE_GREY_BODY: grey-body
         do nothing */
      break;
    case 2:
      /* SURFACE_BC = MO_ATMOSPHERE_TYPE_ZAHNLE: steam atmosphere
         do nothing */
      break;
    case 3:
      /* SURFACE_BC = MO_ATMOSPHERE_TYPE_VOLATILES: self-consistent atmosphere evolution
           with CO2 and H2O volatile using plane-parallel radiative equilibrium model
           of Abe and Matsui (1985)
         do nothing */
      break;
    case 4:
      // MO_ATMOSPHERE_TYPE_HEAT_FLUX: heat flux (prescribed)
      Ap->surface_bc_value /= C->FLUX;
      break;
    case 5:
      /* SURFACE_BC = MO_ATMOSPHERE_TYPE_ENTROPY: entropy
         do nothing */
      break;
    default:
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported SURFACE_BC value %d provided",Ap->SURFACE_BC);
      break;
  }

  Ap->VISCOUS_MANTLE_COOLING_RATE = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-VISCOUS_MANTLE_COOLING_RATE",&Ap->VISCOUS_MANTLE_COOLING_RATE,NULL);CHKERRQ(ierr);

  /* emissivity is constant for SURFACE_BC != MO_ATMOSPHERE_TYPE_VOLATILES */
  Ap->emissivity0 = 1.0; // non-dimensional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-emissivity0",&Ap->emissivity0,NULL);CHKERRQ(ierr);

  Ap->THERMAL_ESCAPE = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-THERMAL_ESCAPE",&Ap->THERMAL_ESCAPE,NULL);CHKERRQ(ierr);

  Ap->CONSTANT_ESCAPE = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-CONSTANT_ESCAPE",&Ap->CONSTANT_ESCAPE,NULL);CHKERRQ(ierr);

  /* Gravitational constant (m^3/kg/s^2) */
  Ap->bigG = 6.67408E-11;
  Ap->bigG *= C->DENSITY;
  Ap->bigG *= PetscPowScalar( C->TIME, 2.0 );

  /* Stefan-Boltzmann constant (W/m^2K^4) */
  Ap->sigma = 5.670367e-08;
  Ap->sigma /= C->SIGMA;

  /* equilibrium temperature of the planet (K) */
  Ap->teqm = 273.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-teqm",&Ap->teqm,NULL);CHKERRQ(ierr);
  Ap->teqm /= C->TEMP;

  Ap->tsurf_poststep_change = -1; // (K) (negative value is OFF)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-tsurf_poststep_change",&Ap->tsurf_poststep_change,NULL);CHKERRQ(ierr);
  Ap->tsurf_poststep_change /= C->TEMP;

  /* for radiative boundary condition at the top surface
     dT = param_utbl_const * [Surface temperature]**3 */
  Ap->PARAM_UTBL = PETSC_TRUE;
  Ap->param_utbl_const = 1.0e-7;
  ierr = PetscOptionsGetBool(NULL,NULL,"-PARAM_UTBL",&Ap->PARAM_UTBL,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-param_utbl_const",&Ap->param_utbl_const,NULL);CHKERRQ(ierr);
  if (Ap->PARAM_UTBL){
      Ap->param_utbl_const *= PetscSqr(C->TEMP);
  } else {
      Ap->param_utbl_const = 0.0;
  }

  /* below here are only used for SURFACE_BC = MO_ATMOSPHERE_TYPE_VOLATILES */

  /* atmosphere reference pressure (Pa) */
  Ap->P0 = 101325.0; // Pa (= 1 atm)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-P0",&Ap->P0,NULL);CHKERRQ(ierr);
  Ap->P0 /= C->PRESSURE;

  /* TODO: check with PS as to whether this is the best approach */
  Ap->gravity_ptr = &P->gravity;
  Ap->radius_ptr = &P->radius;
  Ap->VOLATILE_ptr = &C->VOLATILE;

  /* FIXME the gas constant above is scaled differently for the viscosity laws added by Rob Spaargaren
     this one below is for computing the 1-D atmosphere structure. Should merge together and make consistent */
  Ap->Rgas = 8.3144598; // gas constant (J/K/mol)
  Ap->Rgas *= C->TEMP / C->ENERGY;

  /* Boltzmann constant (J/K) */
  Ap->Avogadro = 6.02214076E23; // 1/mol
  Ap->kB = Ap->Rgas / Ap->Avogadro;

  /* Look for command-line option to determine number of volatiles
     and options prefix for each, e.g -volatile_names CO2,H2O */
  Ap->n_volatiles = 0;
  {
    char      *prefixes[SPIDER_MAX_VOLATILE_SPECIES];
    PetscInt  n_volatiles = SPIDER_MAX_VOLATILE_SPECIES;
    PetscBool set;

    ierr = PetscOptionsGetStringArray(NULL,NULL,"-volatile_names",prefixes,&n_volatiles,&set);CHKERRQ(ierr);
    if (set) {
      PetscInt v;

      Ap->n_volatiles = n_volatiles;
      for (v=0; v<Ap->n_volatiles; ++v) {
        ierr = PetscStrncpy(Ap->volatile_parameters[v].prefix,prefixes[v],sizeof(Ap->volatile_parameters[v].prefix));CHKERRQ(ierr);
        ierr = PetscFree(prefixes[v]);CHKERRQ(ierr);
      }
    }
  }

  /* Get command-line values for all volatiles species */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = VolatileParametersSetFromOptions(&Ap->volatile_parameters[v], C);CHKERRQ(ierr);
  }

  /* Reactions: look for command-line options to determine the number of reactions
     and options for each. These include named reactions and "simple" reactions.
     These are defined with respect to the prefixes for the volatiles, which are translated to integer ids */
  Ap->n_reactions = 0;

  /* Special "named" reactions */
  {
    PetscBool flg;

    ierr = PetscOptionsGetBool(NULL,NULL,"-reaction_ammonia1",NULL,&flg);CHKERRQ(ierr);
    if (flg) {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source",SPIDER_MAX_REACTIONS);
        ierr = ReactionParametersCreateAmmonia1(&Ap->reaction_parameters[Ap->n_reactions],Ap);CHKERRQ(ierr);
        ++Ap->n_reactions;
    }

    ierr = PetscOptionsGetBool(NULL,NULL,"-reaction_carbondioxide1",NULL,&flg);CHKERRQ(ierr);
    if (flg) {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source",SPIDER_MAX_REACTIONS);
        ierr = ReactionParametersCreateCarbonDioxide1(&Ap->reaction_parameters[Ap->n_reactions],Ap);CHKERRQ(ierr);
        ++Ap->n_reactions;
    }

    ierr = PetscOptionsGetBool(NULL,NULL,"-reaction_methane1",NULL,&flg);CHKERRQ(ierr);
    if (flg) {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source",SPIDER_MAX_REACTIONS);
        ierr = ReactionParametersCreateMethane1(&Ap->reaction_parameters[Ap->n_reactions],Ap);CHKERRQ(ierr);
        ++Ap->n_reactions;
    }

    ierr = PetscOptionsGetBool(NULL,NULL,"-reaction_water1",NULL,&flg);CHKERRQ(ierr);
    if (flg) {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source",SPIDER_MAX_REACTIONS);
        ierr = ReactionParametersCreateWater1(&Ap->reaction_parameters[Ap->n_reactions],Ap);CHKERRQ(ierr);
        ++Ap->n_reactions;
    }

    ierr = PetscOptionsGetBool(NULL,NULL,"-reaction_simplewater1",NULL,&flg);CHKERRQ(ierr);
    if (flg) {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source",SPIDER_MAX_REACTIONS);
        ierr = ReactionParametersCreateSimpleWater1(&Ap->reaction_parameters[Ap->n_reactions],Ap);CHKERRQ(ierr);
        ++Ap->n_reactions;
    }

    ierr = PetscOptionsGetBool(NULL,NULL,"-reaction_simplewater2",NULL,&flg);CHKERRQ(ierr);
    if (flg) {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source",SPIDER_MAX_REACTIONS);
        ierr = ReactionParametersCreateSimpleWater2(&Ap->reaction_parameters[Ap->n_reactions],Ap);CHKERRQ(ierr);
        ++Ap->n_reactions;
    }

    ierr = PetscOptionsGetBool(NULL,NULL,"-reaction_simplewater3",NULL,&flg);CHKERRQ(ierr);
    if (flg) {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source",SPIDER_MAX_REACTIONS);
        ierr = ReactionParametersCreateSimpleWater3(&Ap->reaction_parameters[Ap->n_reactions],Ap);CHKERRQ(ierr);
        ++Ap->n_reactions;
    }

  }

  Ap->OXYGEN_FUGACITY = OXYGEN_FUGACITY_NONE;
  {
    PetscInt  OXYGEN_FUGACITY = 0;
    PetscBool OXYGEN_FUGACITYset = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL,NULL,"-OXYGEN_FUGACITY",&OXYGEN_FUGACITY,&OXYGEN_FUGACITYset);CHKERRQ(ierr);
    /* only set OXYGEN_FUGACITY if there are reactions, since fO2 is
       used to also adjust the total atmosphere pressure */
    if( OXYGEN_FUGACITYset && Ap->n_reactions ) Ap->OXYGEN_FUGACITY = OXYGEN_FUGACITY;
  }

  /* radiogenic heating */
  /* aluminium 26 */
  al26->t0 = 0.0; // years
  ierr = PetscOptionsGetScalar(NULL,NULL,"-al26_t0",&al26->t0,NULL);CHKERRQ(ierr);
  al26->t0 /= C->TIMEYRS;
  al26->abundance = 0.0; // fractional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-al26_abundance",&al26->abundance,NULL);CHKERRQ(ierr);
  al26->concentration = 0.0; // ppm
  ierr = PetscOptionsGetScalar(NULL,NULL,"-al_concentration",&al26->concentration,NULL);CHKERRQ(ierr);
  al26->heat_production = 0.3583; // W/kg (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-al26_heat_production",&al26->heat_production,NULL);CHKERRQ(ierr);
  al26->heat_production /= C->HEATGEN;
  al26->half_life = 0.717E6; // years (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-al26_half_life",&al26->half_life,NULL);CHKERRQ(ierr);
  al26->half_life /= C->TIMEYRS;

  /* potassium 40 */
  k40->t0 = 4.55E9; // years
  ierr = PetscOptionsGetScalar(NULL,NULL,"-k40_t0",&k40->t0,NULL);CHKERRQ(ierr);
  k40->t0 /= C->TIMEYRS;
  k40->abundance = 1.1668E-4; // fractional (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-k40_abundance",&k40->abundance,NULL);CHKERRQ(ierr);
  k40->concentration = 0.0; // ppm
  ierr = PetscOptionsGetScalar(NULL,NULL,"-k_concentration",&k40->concentration,NULL);CHKERRQ(ierr);
  k40->heat_production = 2.8761E-5; // W/kg (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-k40_heat_production",&k40->heat_production,NULL);CHKERRQ(ierr);
  k40->heat_production /= C->HEATGEN;
  k40->half_life = 1248E6; // years (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-k40_half_life",&k40->half_life,NULL);CHKERRQ(ierr);
  k40->half_life /= C->TIMEYRS;

  /* iron 60 */
  fe60->t0 = 0.0; // years
  ierr = PetscOptionsGetScalar(NULL,NULL,"-fe60_t0",&fe60->t0,NULL);CHKERRQ(ierr);
  fe60->t0 /= C->TIMEYRS;
  fe60->abundance = 0.0; // fractional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-fe60_abundance",&fe60->abundance,NULL);CHKERRQ(ierr);
  fe60->concentration = 0.0; // ppm
  ierr = PetscOptionsGetScalar(NULL,NULL,"-fe_concentration",&fe60->concentration,NULL);CHKERRQ(ierr);
  fe60->heat_production = 3.6579E-2; // W/kg (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-fe60_heat_production",&fe60->heat_production,NULL);CHKERRQ(ierr);
  fe60->heat_production /= C->HEATGEN;
  fe60->half_life = 2.62E6; // years (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-fe60_half_life",&fe60->half_life,NULL);CHKERRQ(ierr);
  fe60->half_life /= C->TIMEYRS;

  /* thorium 232 */
  th232->t0 = 4.55E9; // years
  ierr = PetscOptionsGetScalar(NULL,NULL,"-th232_t0",&th232->t0,NULL);CHKERRQ(ierr);
  th232->t0 /= C->TIMEYRS;
  th232->abundance = 1.0; // fractional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-th232_abundance",&th232->abundance,NULL);CHKERRQ(ierr);
  th232->concentration = 0.0; // ppm
  ierr = PetscOptionsGetScalar(NULL,NULL,"-th_concentration",&th232->concentration,NULL);CHKERRQ(ierr);
  th232->heat_production = 2.6368E-5; // W/kg (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-th232_heat_production",&th232->heat_production,NULL);CHKERRQ(ierr);
  th232->heat_production /= C->HEATGEN;
  th232->half_life = 14000E6; // years (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-th232_half_life",&th232->half_life,NULL);CHKERRQ(ierr);
  th232->half_life /= C->TIMEYRS;

  /* uranium 235 */
  u235->t0 = 4.55E9; // years
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u235_t0",&u235->t0,NULL);CHKERRQ(ierr);
  u235->t0 /= C->TIMEYRS;
  u235->abundance = 0.0072045; // fractional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u235_abundance",&u235->abundance,NULL);CHKERRQ(ierr);
  u235->concentration = 0.0; // ppm
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u_concentration",&u235->concentration,NULL);CHKERRQ(ierr);
  u235->heat_production = 5.68402E-4; // W/kg (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u235_heat_production",&u235->heat_production,NULL);CHKERRQ(ierr);
  u235->heat_production /= C->HEATGEN;
  u235->half_life = 704E6; // years (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u235_half_life",&u235->half_life,NULL);CHKERRQ(ierr);
  u235->half_life /= C->TIMEYRS;

  /* uranium 238 */
  u238->t0 = 4.55E9; // years
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u238_t0",&u238->t0,NULL);CHKERRQ(ierr);
  u238->t0 /= C->TIMEYRS;
  u238->abundance = 0.0; // fractional
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u238_abundance",&u238->abundance,NULL);CHKERRQ(ierr);
  u238->concentration = 0.0; // ppm
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u_concentration",&u238->concentration,NULL);CHKERRQ(ierr);
  u238->heat_production = 9.4946E-5; // W/kg (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u238_heat_production",&u238->heat_production,NULL);CHKERRQ(ierr);
  u238->heat_production /= C->HEATGEN;
  u238->half_life = 4468E6; // years (Ruedas, 2017)
  ierr = PetscOptionsGetScalar(NULL,NULL,"-u238_half_life",&u238->half_life,NULL);CHKERRQ(ierr);
  u238->half_life /= C->TIMEYRS;

#if 0
  /* compositional parameters */
  P->COMPOSITION = PETSC_FALSE;
  Compp->Brg_initial_fraction = 0.794;
  Compp->Res_Brg_mass_ratio = 1.2362;
  ierr = PetscOptionsGetBool(NULL,NULL,"-COMPOSITION",&P->COMPOSITION,NULL);CHKERRQ(ierr);
  if( P->COMPOSITION ){
    ierr = PetscOptionsGetScalar(NULL,NULL,"-Brg_initial_fraction",&Compp->Brg_initial_fraction,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-Res_Brg_mass_ratio",&Compp->Res_Brg_mass_ratio,NULL);CHKERRQ(ierr);
  }
  Compp->BSE_Brg_mass_ratio_at_liquidus = get_BSE_Brg_mass_ratio( Compp->Brg_initial_fraction, Compp->Res_Brg_mass_ratio);

  /* FIXME below should be initialised elsewhere */
  //Comp->mo_bridgmanite_fraction = -1.0; // updated by code, but initialise here
  //Comp->mo_mass_ratio = -1.0; // updated by code, but initialise here
#endif

  PetscFunctionReturn(0);
}

static PetscErrorCode set_start_time_from_file( Parameters *P , const char * filename )
{

    PetscErrorCode   ierr;
    cJSON            *json=NULL, *time;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_start_time_from_file()\n");CHKERRQ(ierr);

    ierr = read_JSON_file_to_JSON_object( filename, &json );

    /* time from this restart JSON must be passed to the time stepper
       to continue the integration and keep track of absolute time */
    time = cJSON_GetObjectItem(json,"time");
    P->t0 = time->valuedouble;

    cJSON_Delete( json );

    PetscFunctionReturn(0);

}

PetscErrorCode PrintParameters(Parameters const *P)
{
  PetscErrorCode             ierr;
  PetscInt                   i;
  Constants const            *C  = &P->constants;
  AtmosphereParameters const *Ap = &P->atmosphere_parameters;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"                                                                                                    );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"**************** Magma Ocean | Parameters **************\n\n"                                          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15s %s\n"                 ,"[Scaling]"  ,"","Value"                       ,"Units"       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Radius"     ,"",(double)C->RADIUS             ,"m"           );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Temperature","",(double)C->TEMP               ,"K"           );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Entropy"    ,"",(double)C->ENTROPY            ,"J/kg-K"      );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Density"    ,"",(double)C->DENSITY            ,"kg/m^3"      );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s (%.6g years)\n","Time"       ,"",(double)C->TIME               ,"s",(double)C->TIMEYRS);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Area"       ,"",(double)C->AREA               ,"m^2"         );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Volume"     ,"",(double)C->VOLUME             ,"m^3"         );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Mass"       ,"",(double)C->MASS               ,"kg"          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"SEnergy"    ,"",(double)C->SENERGY            ,"J/kg"        );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Energy"     ,"",(double)C->ENERGY             ,"J"           );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Pressure"   ,"",(double)C->PRESSURE           ,"Pa"          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Power"      ,"",(double)C->POWER              ,"W"           );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Flux"       ,"",(double)C->FLUX               ,"W/m^2"       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dP/dR"      ,"",(double)C->DPDR               ,"Pa/m"        );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Gravity"    ,"",(double)C->GRAVITY            ,"m/s^2"       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Kappa"      ,"",(double)C->KAPPA              ,"m^2/s"       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dT/dP"      ,"",(double)C->DTDP               ,"K/Pa"        );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dS/dR"      ,"",(double)C->DSDR               ,"J/kg-K-m"    );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dT/dR"      ,"",(double)C->DTDR               ,"K/m"         );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"GSuper"     ,"",(double)C->GSUPER             ,"K/s^2"       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Visc"       ,"",(double)C->VISC               ,"Pa-s"        );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Log10 Visc" ,"",(double)C->LOG10VISC          ,"log10(Pa-s)" );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Cond"       ,"",(double)C->COND               ,"W/m-K"       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Sigma"      ,"",(double)C->SIGMA              ,"W/m^2-K^4"   );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Lhs"        ,"",(double)C->LHS                ,"kg-K"        );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Rhs"        ,"",(double)C->RHS                ,"W/kg-K"      );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"                                                                                                    );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15s %s\n"    ,"[Parameter]","Non-dim. Value","Dim. Value"                 ,"Units"       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15.6g %-15.6g %s\n","dtmacro"    ,(double)P->dtmacro     ,(double)(P->dtmacro*C->TIME)  ,"s"   );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15d\n"             ,"nstepsmacro",P->nstepsmacro                                               );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15d\n"             ,"numpts_b"   ,P->numpts_b                                                  );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15d\n"             ,"numpts_s"   ,P->numpts_s                                                  );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15.6g %-15.6g %s\n","ic_adiabat_entropy"     ,(double)P->ic_adiabat_entropy       ,(double)(P->ic_adiabat_entropy*C->ENTROPY) ,"J/kg-K"      );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"liquidus data file"         ,P->liquidusFilename                          );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"solidus data file"          ,P->solidusFilename                           );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"alphaSol data file"         ,P->eos2_parameters.lookup.alpha_filename     );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"alphaMel data file"         ,P->eos1_parameters.lookup.alpha_filename     );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"cpSol data file"            ,P->eos2_parameters.lookup.cp_filename        );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"cpMel data file"            ,P->eos1_parameters.lookup.cp_filename        );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"dtdpsSol data file"         ,P->eos2_parameters.lookup.dTdPs_filename     );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"dtdpsMel data file"         ,P->eos1_parameters.lookup.dTdPs_filename     );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"rhoSol data file"           ,P->eos2_parameters.lookup.rho_filename       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"rhoMel data file"           ,P->eos1_parameters.lookup.rho_filename       );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"tempSol data file"          ,P->eos2_parameters.lookup.temp_filename      );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"tempMel data file"          ,P->eos1_parameters.lookup.temp_filename      );CHKERRQ(ierr);
  if (Ap->n_volatiles > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n[Volatile] prefix/name\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                          );CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"%-10D %-15s (.. additional parameters omitted ..)\n",i,Ap->volatile_parameters[i].prefix);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                          );CHKERRQ(ierr);
  }
  if (P->postStepActive) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                          );CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PostStep logic active\n"                                                                             );CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %s\n"                ,"Output Directory"           ,P->outputDirectory                           );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************************\n"                                          );CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode AtmosphereParametersDestroy(AtmosphereParameters* Ap)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;
  for (i=0; i<Ap->n_reactions; ++i) {
    ierr = ReactionParametersDestroy(&Ap->reaction_parameters[i]);CHKERRQ(ierr);
  }
  Ap->n_reactions = 0;
  PetscFunctionReturn(0);
}

PetscErrorCode ParametersDestroy(Parameters* parameters)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = AtmosphereParametersDestroy(&parameters->atmosphere_parameters);CHKERRQ(ierr);
  EosParametersDestroy(parameters);
  PetscFunctionReturn(0);
}
