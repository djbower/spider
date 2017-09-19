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
    // 35 constants to set
    SQRTST = PetscSqrtScalar( ENTROPY * TEMPERATURE );

    C->RADIUS   = RADIUS;
    C->TEMP     = TEMPERATURE;
    C->ENTROPY  = ENTROPY;
    C->DENSITY  = DENSITY;
    C->AREA     = PetscSqr( C->RADIUS );
    C->AREAG    = C->AREA * 4.0 * PETSC_PI;
    C->VOLUME   = C->AREA * C->RADIUS;
    C->VOLUMEG  = C->VOLUME * 4.0 * PETSC_PI;
    C->MASS     = C->DENSITY * C->VOLUME;
    C->MASSG    = C->MASS * 4.0 * PETSC_PI;
    C->TIME     = C->RADIUS / SQRTST;
    C->TIMEYRS  = C->TIME / (60.0*60.0*24.0*365.25);
    C->SENERGY  = C->ENTROPY * C->TEMP;
    C->ENERGY   = C->SENERGY * C->MASS;
    C->ENERGYG  = C->ENERGY * 4.0 * PETSC_PI;
    C->PRESSURE = C->ENTROPY * C->TEMP * C->DENSITY;
    C->POWER    = C->ENERGY / C->TIME;
    C->POWERG   = C->POWER * 4.0 * PETSC_PI;
    C->FLUX     = C->POWER / C->AREA;
    C->DPDR     = C->PRESSURE / C->RADIUS;
    C->ALPHA    = 1.0 / C->TEMP;
    C->GRAVITY  = (C->ENTROPY * C->TEMP) / C->RADIUS;
    C->KAPPA    = C->RADIUS * SQRTST;
    C->DTDP     = 1.0 / (C->DENSITY * C->ENTROPY);
    C->DSDR     = C->ENTROPY / C->RADIUS;
    C->DTDR     = C->TEMP / C->RADIUS;
    C->GSUPER   = C->GRAVITY * C->DTDR;
    C->ETA      = C->DENSITY * C->KAPPA;
    C->LOG10ETA = PetscLog10Real( C->ETA );
    C->NU       = C->KAPPA;
    C->COND     = C->ENTROPY * C->DENSITY * C->KAPPA;
    C->SIGMA    = C->FLUX * 1.0 / PetscPowScalar( C->TEMP, 4.0 );
    C->LHS      = C->DENSITY * C->VOLUME * C->TEMP;
    C->LHSG     = C->LHS * 4.0 * PETSC_PI;
    C->RHS      = C->ENTROPY / C->TIME;

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
    ierr = SetConstants(C,RADIUS0,TEMPERATURE0,ENTROPY0,DENSITY0);CHKERRQ(ierr); /* see global_defs.h */
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
  Constants const      *C  = &P->constants;

  PetscFunctionBegin;
  ierr = InitializeConstantsAndSetFromOptions(&P->constants);CHKERRQ(ierr);

  //TODO these should really happen at the same time, and be included here.

  /* Default discretization parameters */
  P->nstepsmacro = 1000;
  P->maxsteps    = 100000000; /* Effectively infinite */
  P->dtmacro     = 1000000.0;
  P->t0          = 0.0;

  P->numpts_b = 200;
//P->numpts_b= 278
//P->numpts_b= 372
//P->numpts_b= 656
//P->numpts_b= 2939
//P->numpts_b= 5802
//P->numpts_b= 11532
//P->numpts_b= 22996
//P->numpts_b= 45928
  P->numpts_s = P->numpts_b + 1;

  P->monitor = PETSC_TRUE;

  /* For each entry in parameters, we set a default value and immediately scale it.
     Dimensional/unscaled quantities are not explicitly stored. */
  // TODO

  // TODO -- pasted begin
    // all SI units unless non-dimensional

    // initial entropy at top of adiabat (J/kgK)
    P->sinit = 3052.885602072091;
    ierr = PetscOptionsGetScalar(NULL,NULL,"-sinit",&P->sinit,NULL);CHKERRQ(ierr);
    P->sinit /= C->ENTROPY;

    // initial entropy gradient (J/kgKm)
    P->ic_dsdr = -4.6978890285209187e-07;
    P->ic_dsdr /= C->DSDR;
    // radius of planet (m)
    P->radius = 6371000.0;
    P->radius /= C->RADIUS;
    // core size (non-dimenisonal)
    P->coresize = 0.55;
    // surface density (kg/m^3) for Adams-Williamson EOS for pressure
    P->rhos = 4078.95095544;
    P->rhos /= C->DENSITY;
    // parameter (1/m) for Adams-Williamson EOS for pressure
    P->beta = 1.1115348931000002e-07;
    P->beta *= C->RADIUS;
    // grain size (m)
    P->grain = 1.0E-3;
    P->grain /= C->RADIUS;
    // gravity (m/s^2), must be negative
    P->gravity = -10.0;
    P->gravity /= C->GRAVITY;
    // melt fraction threshold for rheology
    P->phi_critical = 0.4; // non dimensional
    // melt fraction transition width for rheology
    /* when PHI_WIDTH = 0.15, oscillations appear in the volatile contents
   due to a rigid crust forming at the top of the model.  Reducing the value
   to 0.2 helps to alleviate this problem.  So evidently the viscosity contrast
   across nodes matters. */
    P->phi_width = 0.15; // non dimensional
    // melt fraction shape transition for skew
    P->phi_skew = 0.0; // non dimensional
    // core density (kg/m^3)
    P->rho_core = 10738.332568062382;
    P->rho_core /= C->DENSITY;
    // heat capacity of core (J/kgK)
    P->cp_core = 880.0;
    P->cp_core /= C->ENTROPY;
    // mass-weighted average core temperature as a fraction
    // of CMB temperature (non-dimensional)
    P->tfac_core_avg = 1.147;
    // smoothing width
    P->swidth = 1.0E-2;
    // solid viscosity (Pa.s)
    P->log10visc_sol = 21.0;
    P->log10visc_sol -= C->LOG10ETA;
    // solid conductivity (W/mK)
    P->cond_sol = 4.0;
    P->cond_sol /= C->COND;
    // melt viscosity (Pa.s)
    P->log10visc_mel = 2.0;
    P->log10visc_mel -= C->LOG10ETA;
    // melt conductivity (W/mK)
    P->cond_mel = 4.0;
    P->cond_mel /= C->COND;

    /* atmosphere parameters
         MODEL = MO_ATMOSPHERE_TYPE_GREY_BODY: grey-body
         MODEL = MO_ATMOSPHERE_TYPE_ZAHNLE: FIXME: currently broken
         MODEL = MO_ATMOSPHERE_TYPE_VOLATILES: self-consistent atmosphere evolution
                  with CO2 and H2O volatiles
                  uses plane-parallel radiative eqm model
                  of Abe and Matsui (1985)
    */
    Ap->MODEL=MO_ATMOSPHERE_TYPE_GREY_BODY;
    /* for legacy purposes
       if HYBRID is set, then the boundary condition will switch
       to the upper mantle cooling rate once the rheological
       transition is reached.  This prevents a lid from forming at
       the top of the model.
       TODO: this implies that the emissivity is around 1.0E-7,
       which is unphysical unless you are appealing to a massive
       massive atmosphere */
    Ap->HYBRID=0;
    // emissivity is constant for MODEL != MO_ATMOSPHERE_TYPE_VOLATILES
    Ap->EMISSIVITY0 = 1.0;
    // Stefan-Boltzmann constant (W/m^2K^4)
    Ap->SIGMA = 5.670367e-08;
    // equilibrium temperature of the planet (K)
    Ap->TEQM = 273.0;
    /* for radiative boundary condition at the top surface
       dT = constbc * [Surface temperature]**3
       FIXME: need to non-dimensionalise below! */
    //Ap->CONSTBC = 1.0e-07; // FIXME check units
    Ap->CONSTBC = 0;

    /* below here are only used for MODEL = MO_ATMOSPHERE_TYPE_VOLATILES */
    // FIXME: for convenience at the moment, duplicate these values */
    Ap->RADIUS = P->radius; // FIXME: dimensional
    Ap->GRAVITY = P->gravity; // FIXME: dimensional
    /* VOLSCALE enables us to scale the volatile equations to the same
       order of magnitude as entropy, and thus ensure that the residual
       based on the solution vector is not biased. */
    // FIXME: need to figure out scalings of this ODE
    //Ap->VOLSCALE = 1.0E2; // wt %
    Ap->VOLSCALE = 1.0E6; // ppm
    // initial volatile contents in the liquid magma ocean
    Ap->H2O_INITIAL = 0.0;
    Ap->CO2_INITIAL = 0.0;
    // Elkins-Tanton (2008) case 1
    //P->h2o_initial = 500.0;
    //P->co2_initial = 100.0;
    // Elkins-Tanton (2008) case 2
    //P->h2o_initial = 5000.0;
    //P->co2_initial = 1000.0;
    // Elkins-Tanton (2008) case 3
    //P->h2o_initial = 0.0;
    //P->co2_initial = 6000.0;
    // atmosphere reference pressure (Pa)
    Ap->P0 = 101325.0; // Pa (= 1 atm)
    // distribution coefficients are given in supplementary material of ET08
    // distribution coefficient between solid and melt (non-dimensional)
    Ap->H2O_KDIST = 1.0E-4;
    // TODO: water saturation limit of 10 ppm?
    // absorption (m^2/kg)
    Ap->H2O_KABS = 0.01;
    // next two from Lebrun et al. (2013)
    Ap->H2O_HENRY = 6.8E-8; // must be mass fraction/Pa
    Ap->H2O_HENRY_POW = 1.4285714285714286; // (1.0/0.7)
    // distribution coefficient between solid and melt (non-dimensional)
    Ap->CO2_KDIST = 5.0E-4;
    // TODO: water saturation limit of 0.03 ppm
    // absorption (m^2/kg)
    Ap->CO2_KABS = 0.05;
    // next two from Lebrun et al. (2013)
    Ap->CO2_HENRY = 4.4E-12; // must be mass fraction/Pa
    Ap->CO2_HENRY_POW = 1.0;

  // TODO -- pasted end

    /* TODO: non-dimensionalise other atmosphere parameters */
    /* SIGMA and TEQM should get us up and running for the standard
       grey-body case */
    Ap->SIGMA /= C->SIGMA;
    Ap->TEQM /= C->TEMP;
    Ap->CONSTBC /= 1.0; // FIXME to UPDATE!

  /* Additional Parameters for Atmosphere */
  // TODO

    // TODO move this up so it's all done at the same time..
  /* Get lookup tables */
  ierr = SetLookups(P);CHKERRQ(ierr);

  /* Get grid parameters */
  {
    PetscBool set = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&P->numpts_s,&set);CHKERRQ(ierr);
    if (set) P->numpts_b = P->numpts_s + 1;
  }

  /* Obtain command-line options for simulation time frame and monitoring */
  {
    PetscReal dtmacro_years;
    PetscBool dtmacro_set = PETSC_FALSE, dtmacro_years_set = PETSC_FALSE, nstepsmacro_set = PETSC_FALSE,
              early = PETSC_FALSE, middle=PETSC_FALSE, late = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,NULL,"-monitor",&P->monitor,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-nstepsmacro",&P->nstepsmacro,&nstepsmacro_set);CHKERRQ(ierr);
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
      if (early) {
        /* early evolution to about 10 kyr */
        P->nstepsmacro = 1000;
        P->dtmacro = 1000000;
      } else if (middle) {
        /* middle evolution to about 100 Myr */
        P->nstepsmacro = 10000;
        P->dtmacro = 1000000000;
      } else if (late) {
        /* late evolution to 4.55 Byr */
        P->nstepsmacro = 455;
        P->dtmacro = 1000000000000;
      }
    } else {
      if (dtmacro_set && dtmacro_years_set) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: both -dtmacro and -dtmacro_years provided. Using -dtmacro\n");CHKERRQ(ierr);
      } else if (dtmacro_years_set) {
        P->dtmacro = dtmacro_years / C->TIMEYRS;
      }
    }
  }

  /* Additional Parameters for Atmosphere */
  // TODO

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
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Area_G"     ,"",C->AREAG,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Volume"     ,"",C->VOLUME,"m^3");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Volume_G"   ,"",C->VOLUMEG,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Mass"       ,"",C->MASS,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Mass_G"     ,"",C->MASSG,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"SEnergy"    ,"",C->SENERGY,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Energy"     ,"",C->ENERGY,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Energy_G"   ,"",C->ENERGYG,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Pressure"   ,"",C->PRESSURE,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Power"      ,"",C->POWER,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Power_G"    ,"",C->POWERG,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Flux"       ,"",C->FLUX,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dPdR"       ,"",C->DPDR,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Alpha"      ,"",C->ALPHA,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Gravity"    ,"",C->GRAVITY,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Kappa"      ,"",C->KAPPA,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dTdP"       ,"",C->DTDP,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dSdR"       ,"",C->DSDR,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"dTdR"       ,"",C->DTDR,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"GSuper"     ,"",C->GSUPER,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Eta"        ,"",C->ETA,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Log10 Eta"  ,"",C->LOG10ETA,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"NU"         ,"",C->NU,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Cond"       ,"",C->COND,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Sigma"      ,"",C->SIGMA,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Lhs"        ,"",C->LHS,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Lhs_G"      ,"",C->LHSG,"");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15.6g %-6s\n"             ,"Rhs"        ,"",C->RHS,"");CHKERRQ(ierr);
  // TODO finish adding units, do other cleanup

  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"                                                                                    );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15s %-15s %s\n"    ,"[Parameter]","Non-dim. Value","Dim. Value"       ,"Units" );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15.6g %-15.6g %s\n","dtmacro"    ,P->dtmacro      ,P->dtmacro*C->TIME ,"s"     );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15d\n"             ,"nstepsmacro",P->nstepsmacro                               );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15d\n"             ,"numpts_b"   ,P->numpts_b                                  );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-15s %-15.6g %-15.6g %s\n","S_init"     ,P->sinit        ,P->sinit*C->ENTROPY,"J/kg-K");CHKERRQ(ierr);
  // TODO add the rest..
  ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                            );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************************\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
