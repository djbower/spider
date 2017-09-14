/*
Parameter Management

Parameters should only ever be set by the functions in this file. That is, everywhere else they should be considered read-only.

Custom PETSc command line options should only ever be parsed here.

 */

#include "parameters.h"
#include "ctx.h"

PetscErrorCode set_constants( Constants *C )
{
    PetscScalar SQRTST;

    PetscFunctionBeginUser;

    // 35 constants to set
    SQRTST = PetscSqrtScalar( ENTROPY0 * TEMPERATURE0 );

    C->RADIUS   = RADIUS0;
    C->TEMP     = TEMPERATURE0;
    C->ENTROPY  = ENTROPY0;
    C->DENSITY  = DENSITY0;
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
// TODO - paste end ...

PetscErrorCode InitializeParameters(Parameters *P) 
{
  AtmosphereParameters *Ap = &P->atmosphere_parameters;
  Constants            *C  = &P->constants;

  PetscFunctionBeginUser;

  /* Set default scalings */
  // TODO
  // TODO paste begin ..
  set_constants(&P->constants);
  // TODO paste end ..

  /* For each entry in parameters, we set a default value and immediately scale it. Dimensional/unscaled quantities are not explicitly stored. */
  // TODO

  // TODO -- pasted begin
    // all SI units unless non-dimensional

    // initial entropy at top of adiabat (J/kgK)
    P->sinit = 3052.885602072091;
    // initial entropy gradient (J/kgKm)
    P->ic_dsdr = -4.6978890285209187e-07;
    // radius of planet (m)
    P->radius = 6371000.0;
    // core size (non-dimenisonal)
    P->coresize = 0.55;
    // surface density (kg/m^3) for Adams-Williamson EOS for pressure
    P->rhos = 4078.95095544;
    // parameter (1/m) for Adams-Williamson EOS for pressure
    P->beta = 1.1115348931000002e-07;
    // grain size (m)
    P->grain = 1.0E-3;
    // gravity (m/s^2), must be negative
    P->gravity = -10.0;
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
    // heat capacity of core (J/kgK)
    P->cp_core = 880.0;
    // mass-weighted average core temperature as a fraction
    // of CMB temperature (non-dimensional)
    P->tfac_core_avg = 1.147;
    // smoothing width
    P->swidth = 1.0E-2;
    // solid viscosity (Pa.s)
    P->log10visc_sol = 21.0;
    // solid conductivity (W/mK)
    P->cond_sol = 4.0;
    // melt viscosity (Pa.s)
    P->log10visc_mel = 2.0;
    // melt conductivity (W/mK)
    P->cond_mel = 4.0;

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

    // TODO - pasted these in..
    // We want to have in-place application of these..
    P->sinit /= C->ENTROPY;
    P->ic_dsdr /= C->DSDR;
    P->radius /= C->RADIUS;
    P->rhos /= C->DENSITY;
    P->beta *= C->RADIUS;
    P->grain /= C->RADIUS;
    P->gravity /= C->GRAVITY;
    P->rho_core /= C->DENSITY;
    P->cp_core /= C->ENTROPY;
    P->log10visc_sol -= C->LOG10ETA;
    P->cond_sol /= C->COND;
    P->log10visc_mel -= C->LOG10ETA;
    P->cond_mel /= C->COND;

    /* TODO: non-dimensionalise other atmosphere parameters */
    /* SIGMA and TEQM should get us up and running for the standard
       grey-body case */
    Ap->SIGMA /= C->SIGMA;
    Ap->TEQM /= C->TEMP;
    Ap->CONSTBC /= 1.0; // FIXME to UPDATE!

  /* Additional Parameters for Atmosphere */
  // TODO

  PetscFunctionReturn(0);
}

PetscErrorCode SetParametersFromOptions(Parameters *parameters)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* Store previous scalings */
  // TODO

  /*  Update Scalings */
  // TODO

  /* For each entry in parameters, we set a default value and immediately re-scale it */
  // TODO

  /* Additional Parameters for Atmosphere */
  // TODO

  PetscFunctionReturn(0);
}

PetscErrorCode PrintParameters(Parameters const *P,FILE *file)
{
  PetscErrorCode ierr;
  Constants const *C = &P->constants;

  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD," ******** Magma Ocean | Parameters ********\n");CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Scalings:\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Radius:   %f m (%f km)\n",C->RADIUS,C->RADIUS/1000.0);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Time:     %f s (%f years)\n",C->TIME,C->TIMEYRS);CHKERRQ(ierr);
  // TODO add the rest..
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Parameter     Non-dim. Value   Dim Value [Units]\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%10s  %12.6f  %12.6f [??]\n","S_init",P->sinit,P->sinit*P->constants.ENTROPY);CHKERRQ(ierr);
  // TODO add the rest..


  ierr = PetscPrintf(PETSC_COMM_WORLD," ******************************************\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
