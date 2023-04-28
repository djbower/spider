/*
Parameter Management

Parameters should only ever be set by the functions in this file. That is, everywhere else they should be considered read-only.

Custom PETSc command line options should only ever be parsed during the populate of the Parameters struct here.
*/

#include "parameters.h"
#include "ctx.h"
#include "eos.h"
#include "eos_composite.h"
#include "interp.h"
#include "util.h"

static PetscErrorCode VolatileParametersCreate(VolatileParameters *);
static PetscErrorCode RadionuclideParametersCreate(RadionuclideParameters *);
static PetscErrorCode AtmosphereParametersSetFromOptions(Parameters, const ScalingConstants, const FundamentalConstants);

static PetscErrorCode ScalingConstantsSet(ScalingConstants SC, PetscReal TEMPERATURE, PetscReal RADIUS, PetscReal TIME, PetscReal PRESSURE, PetscReal VOLATILE)
{
  /* constants used to scale the physical problem are largely chosen based on numerical considerations.
     Factors of 4 pi associated with spherical geometry are excluded, but are reintroduced in output
     routines to give the correct (meaningful) physical values */

  PetscScalar SQRTST;

  PetscFunctionBeginUser;

  /* these 5 scaling constants can be set by the user, and the others
     subsequently derived */
  /* the first four are used to non-dimensionalise the problem, but
     should be chosen to scale solution quantities and their time
     derivatives to around unit */
  SC->TEMP = TEMPERATURE; // K
  SC->RADIUS = RADIUS;    // m
  SC->TIME = TIME;        // s
  /* pressure is only relevant if solving for the atmosphere, and
     is usually around 1-100 bar */
  SC->PRESSURE = PRESSURE; // Pa
  /* volatile scales the volatile mass balance, since mass is already
     defined by the four quantities above we must specify an extra
     scaling to scale the numerical problem appropriately */
  SC->VOLATILE = VOLATILE;
  /* below are derived from above */
  /* note: factors of 4 pi are excluded */
  SC->AREA = PetscSqr(SC->RADIUS);                          // m^2
  SC->VOLUME = SC->AREA * SC->RADIUS;                       // m^3
  SC->ENTROPY = PetscSqr(SC->RADIUS / SC->TIME) / SC->TEMP; // J/kg/K
  SQRTST = PetscSqrtScalar(SC->ENTROPY * SC->TEMP);
  SC->MASS = SC->PRESSURE * SC->VOLUME / (SC->ENTROPY * SC->TEMP); // kg
  SC->DENSITY = SC->PRESSURE / (SC->ENTROPY * SC->TEMP);           // kg/m^3
  SC->TIMEYRS = SC->TIME / (60.0 * 60.0 * 24.0 * 365.25);          // years
  SC->SENERGY = SC->ENTROPY * SC->TEMP;                            // J/kg
  SC->ENERGY = SC->SENERGY * SC->MASS;                             // J
  SC->POWER = SC->ENERGY / SC->TIME;                               // W
  SC->FLUX = SC->POWER / SC->AREA;                                 // W/m^2
  SC->DPDR = SC->PRESSURE / SC->RADIUS;                            // Pa/m
  SC->GRAVITY = SC->ENTROPY * SC->TEMP / SC->RADIUS;               // m/s^2
  SC->KAPPA = SC->RADIUS * SQRTST;                                 // m^2/s
  SC->DSDP = SC->ENTROPY / SC->PRESSURE;                           // K/Pa
  SC->DSDR = SC->ENTROPY / SC->RADIUS;                             // J/kg/K/m
  SC->DTDP = SC->TEMP / SC->PRESSURE;                              // K/Pa
  SC->DTDR = SC->TEMP / SC->RADIUS;                                // K/m
  SC->GSUPER = SC->GRAVITY * SC->DTDR;
  SC->VISC = SC->DENSITY * SC->KAPPA;                                           // Pa.s
  SC->LOG10VISC = PetscLog10Real(SC->VISC);                                     // log10(Pa.s)
  SC->COND = SC->ENTROPY * SC->DENSITY * SC->KAPPA;                             // W/m/K
  SC->SIGMA = SC->FLUX * 1.0 / PetscPowScalar(SC->TEMP, 4.0);                   // W/m^2/K^4
  SC->HEATGEN = PetscPowScalar(SC->ENTROPY * SC->TEMP, 3.0 / 2.0) / SC->RADIUS; // W/kg

  PetscFunctionReturn(0);
}

static PetscErrorCode ScalingConstantsSetFromOptions(ScalingConstants SC)
{
  PetscErrorCode ierr;
  PetscScalar TEMPERATURE0, RADIUS0, TIME0, PRESSURE0, VOLATILE0;

  PetscFunctionBeginUser;

  ierr = PetscOptionsGetPositiveScalar("-temperature0", &TEMPERATURE0, 1.0E3, NULL);
  CHKERRQ(ierr); // K
  ierr = PetscOptionsGetPositiveScalar("-radius0", &RADIUS0, 1.0E6, NULL);
  CHKERRQ(ierr); // m
  ierr = PetscOptionsGetPositiveScalar("-time0", &TIME0, 3.154E7, NULL);
  CHKERRQ(ierr); // s
  ierr = PetscOptionsGetPositiveScalar("-pressure0", &PRESSURE0, 1.0E7, NULL);
  CHKERRQ(ierr); // Pa
  /* volatile0 can be set a priori, since it scales the abundance by mass of a volatile
     to its scaled partial pressure.  Since most Henry coefficients are around 1E-6 to 1E-1
     this value can usually be assumed to be around 1.0E-10 when acccounting for the conversion
     from mass fraction to ppmw which introduces an extra factor of 1E6. */
  ierr = PetscOptionsGetPositiveScalar("-volatile0", &VOLATILE0, 1.0E-10, NULL);
  CHKERRQ(ierr);
  ierr = ScalingConstantsSet(SC, TEMPERATURE0, RADIUS0, TIME0, PRESSURE0, VOLATILE0);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode RadionuclideParametersSetFromOptions(RadionuclideParameters Rp, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  char buf[1024]; /* max size */

  PetscFunctionBeginUser;

  /* Accept -prefix_YYY to populate vp->YYY. Most are required and an error is thrown
     if they are missing. Note that this code has a lot of duplication */
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", Rp->prefix, "_t0");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &Rp->t0, 0.0, NULL);
  CHKERRQ(ierr); // years
  Rp->t0 /= SC->TIMEYRS;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", Rp->prefix, "_abundance");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &Rp->abundance, 0.0, NULL);
  CHKERRQ(ierr); // fractional

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", Rp->prefix, "_concentration");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &Rp->concentration, 0.0, NULL);
  CHKERRQ(ierr);               // ppmw
  Rp->concentration *= 1.0E-6; // to mass fraction

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", Rp->prefix, "_heat_production");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &Rp->heat_production, 0.0, NULL);
  CHKERRQ(ierr); // W/kg
  Rp->heat_production /= SC->HEATGEN;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", Rp->prefix, "_half_life");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &Rp->half_life, 0.0, NULL);
  CHKERRQ(ierr); // years /* undefined problem with zero? */
  Rp->half_life /= SC->TIMEYRS;

  PetscFunctionReturn(0);
}

static PetscErrorCode VolatileParametersSetFromOptions(VolatileParameters vp, const AtmosphereParameters Ap, const ScalingConstants SC, const FundamentalConstants FC)
{
  PetscErrorCode ierr;
  char buf[1024]; /* max size */
  PetscBool set;

  PetscFunctionBeginUser;
  /* Accept -prefix_YYY to populate vp->YYY. Most are required and an error is thrown
     if they are missing. Note that this code has a lot of duplication */
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_initial_total_abundance");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->initial_total_abundance, 0.0, NULL);
  CHKERRQ(ierr);
  vp->initial_total_abundance /= 1.0E6 * SC->VOLATILE;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_initial_ocean_moles");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->initial_ocean_moles, 0.0, NULL);
  CHKERRQ(ierr);

  /* used either to set the pressure to this value for the ic, or alternatively
     used as the initial guess to solve for the pressure to match the initial
     total abundance accounting for reactions */
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_initial_atmos_pressure");
  CHKERRQ(ierr);
  /* default value below is to recover original behaviour, but can be updated */
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->initial_atmos_pressure, 0.1 * SC->PRESSURE, NULL);
  CHKERRQ(ierr);
  vp->initial_atmos_pressure /= SC->PRESSURE;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_kdist");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->kdist, 0.0, NULL);
  CHKERRQ(ierr);

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_kabs");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->kabs, 0.0, NULL);
  CHKERRQ(ierr);
  vp->kabs *= SC->DENSITY * SC->RADIUS;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_SOLUBILITY");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveInt(buf, &vp->SOLUBILITY, 1, NULL);
  CHKERRQ(ierr); // Modified Henry's law

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_henry_pow");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->henry_pow, 1.0, NULL);
  CHKERRQ(ierr);

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_henry");
  CHKERRQ(ierr);
  /* default is no dissolved volatile content */
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->henry, 0.0, NULL);
  CHKERRQ(ierr);
  /* ensure that V->dxdp calculation does not return NaN when vp->henry is zero */
  if (vp->henry == 0)
  {
    vp->henry_pow = 1.0;
  }
  vp->henry /= 1.0E6 * SC->VOLATILE * PetscPowScalar(SC->PRESSURE, -1.0 / vp->henry_pow);

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_molar_mass");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->molar_mass, 0.0, &set);
  CHKERRQ(ierr);
  /* there is no way to pick a suitable default, so force the user to specify */
  if (!set)
    SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "Missing argument: %s", buf);
  /* in principle, we could fully non-dimensionalise using Avogadro's constant, but then
     the per-particle value would be 23 orders of magnitude less than this value */
  vp->molar_mass /= SC->MASS;

  /* escape related */
  vp->jeans_value = 0.0;
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_jeans_value");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &vp->jeans_value, NULL);
  CHKERRQ(ierr);

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_R_thermal_escape_value");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->R_thermal_escape_value, 0.0, NULL);
  CHKERRQ(ierr);

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_constant_escape_value");
  CHKERRQ(ierr);
  /* in principle, could have negative escape to add atmospheric material (e.g. from an external source) */
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->constant_escape_value, 0.0, NULL);
  CHKERRQ(ierr);
  vp->constant_escape_value *= SC->TIME / (SC->VOLATILE * SC->MASS);
  /* recall that the code ignores the 4.0*pi prefactor, and it is reintroduced for output only */
  vp->constant_escape_value /= 4.0 * PETSC_PI;

  /* for Zahnle et al. (2019), Eqn. 3, can precompute many factors to then multiply by
     the volume mixing ratio (of H2) which evolves during the evolution */
  /* get prefactor */
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_R_zahnle_escape_value");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->R_zahnle_escape_value, 0.0, NULL);
  CHKERRQ(ierr);
  /* the value provided by the user is a non-dimensional multipler (must be
     unity to recover the original expression, and zero will turn it off for this
     volatile) */
  /* non-dimensional scaling */
  vp->R_zahnle_escape_value *= PetscSqr(SC->RADIUS) * SC->TIME / (SC->VOLATILE * FC->AVOGADRO);
  /* Zahnle et al. (2019), factors from Eqn 3 */
  /* if solar_xuv_factor eventually time dependent, will need to move calculation
     into the time-stepper */
  vp->R_zahnle_escape_value *= Ap->solar_xuv_factor;
  vp->R_zahnle_escape_value /= PetscSqrtScalar(1.0 + 0.006 * PetscSqr(Ap->solar_xuv_factor));
  vp->R_zahnle_escape_value *= PetscSqr(*Ap->radius_ptr); /* at planetary radius (approx) */
  /* 10^16 since need units of m^-2 and not cm^-2 */
  vp->R_zahnle_escape_value *= PetscPowScalar(10.0, 16) * vp->molar_mass; /* mass flux */
  /* to compute escape flux, just need to multiply above by mixing ratio (of H2),
     which evolves during the evolution */

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_cross_section");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf, &vp->cross_section, 1.0E-18, NULL);
  CHKERRQ(ierr); // m^2, Johnson et al. (2015), N2+N2 collisions
  vp->cross_section /= SC->AREA;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_poststep_change");
  CHKERRQ(ierr);
  vp->poststep_change = -1; // fractional (negative value is OFF)
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &vp->poststep_change, &set);
  CHKERRQ(ierr);

  /* if volatile is a pseudo-volatile, we just prescribe the evolution of its pressure
     as a function of surface temperature, using column data from a file */
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", vp->prefix, "_TP_filename");
  CHKERRQ(ierr);
  ierr = PetscStrcpy(vp->TP_filename, "TP_file.json");
  CHKERRQ(ierr);
  if (Ap->PSEUDO_VOLATILES)
  {
    /* get filename */
    ierr = PetscOptionsGetString(NULL, NULL, buf, vp->TP_filename, PETSC_MAX_PATH_LEN, NULL);
    CHKERRQ(ierr);
    /* create interp1d */
    /* note: second column is log10 of pressure in Pa, so do not scale by SC->PRESSURE */
    ierr = Interp1dCreateAndSet(vp->TP_filename, &vp->TP_interp, SC->TEMP, 1.0);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
This function (and subfunctions) should be the only place that
custom command-line parameters (those not defined by PETSc) should be accessed.

Parameters are specified by the user in dimensional (unscaled) form,
but they are all stored in non-dimensional (scaled) form.
*/

PetscErrorCode ParametersSetFromOptions(Parameters P)
{
  PetscErrorCode ierr;
  FundamentalConstants const FC = P->fundamental_constants;
  ScalingConstants const SC = P->scaling_constants;

  PetscFunctionBegin;

  /* Constants (scalings) must be set first, as they are used to scale
     other parameters */
  /* since this sets constants, cannot use C shorthand above which is read only */
  ierr = ScalingConstantsSetFromOptions(P->scaling_constants);
  CHKERRQ(ierr);

  ierr = FundamentalConstantsSet(P->fundamental_constants, SC);
  CHKERRQ(ierr);

  /* Must set EOS after setting constants, but before boundary conditions
     since EOS might be required to map temperature to entropy
     and vice versa for boundary conditions */

  /* For each entry in parameters, we set a default value and immediately scale it.
     Dimensional/unscaled quantities are not explicitly stored.
     All SI units unless non-dimensional.
     */

  /* Time frame parameters */
  P->maxsteps = 100000000; /* Effectively infinite */

  P->nstepsmacro = 2;
  ierr = PetscOptionsGetInt(NULL, NULL, "-nstepsmacro", &P->nstepsmacro, NULL);
  CHKERRQ(ierr);

  /* initial macro step (starting) */
  /* can also be over-written later if restarting from file */
  P->stepmacro = 0;
  ierr = PetscOptionsGetInt(NULL, NULL, "-stepmacro", &P->stepmacro, NULL);
  CHKERRQ(ierr);

  /* start time (years) P->t0 is set further down, since it may
     acquire a value from a restart value */

  /* step time (years) */
  P->dtmacro = 100;
  ierr = PetscOptionsGetReal(NULL, NULL, "-dtmacro", &P->dtmacro, NULL);
  CHKERRQ(ierr);
  P->dtmacro /= SC->TIMEYRS; // non-dimensional for time stepping

  /* Grid parameters */
  ierr = PetscOptionsGetPositiveInt("-n", &P->numpts_b, 200, NULL);
  CHKERRQ(ierr);
  if (P->numpts_b < 2)
  {
    SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "numpts_b must >= 2 (currently %d)", P->numpts_b);
  }
  P->numpts_s = P->numpts_b - 1;

  /* this allows backward compatibility with older versions of the code
     that worked directly with radius coordinates, rather than mass coordinates */
  P->MASS_COORDINATES = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-MASS_COORDINATES", &P->MASS_COORDINATES, NULL);
  CHKERRQ(ierr);

  /* RollBack and PostStep options */
  P->rollBackActive = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-activate_rollback", &P->rollBackActive, NULL);
  CHKERRQ(ierr);
  P->postStepActive = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-activate_poststep", &P->postStepActive, NULL);
  CHKERRQ(ierr);

  /* Output Options */
  P->monitor = PETSC_TRUE;
  ierr = PetscStrcpy(P->outputDirectory, "output");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, "-outputDirectory", P->outputDirectory, PETSC_MAX_PATH_LEN, NULL);
  CHKERRQ(ierr);

  /* radius of planet (m) */
  ierr = PetscOptionsGetPositiveScalar("-radius", &P->radius, 6371000.0, NULL);
  CHKERRQ(ierr); // m
  P->radius /= SC->RADIUS;

  /* core radius relative to physical radius i.e. radius */
  /* therefore, scaled (code) core radius is P->coresize * P->radius
     and actual physical core radius is P->coresize * P->radius * SC->RADIUS */
  ierr = PetscOptionsGetPositiveScalar("-coresize", &P->coresize, 0.55, NULL);
  CHKERRQ(ierr); // Earth core radius

  ierr = PetscOptionsGetPositiveInt("-mixing_length", &P->mixing_length, 1, NULL);
  CHKERRQ(ierr);
  if (P->mixing_length == 1)
  {
    /* See fig. 2 in Kamata (2018), JGR */
    /* defaults below are for conventional mixing length (distance to nearest boundary) */
    ierr = PetscOptionsGetPositiveScalar("-mixing_length_peak_location", &P->mixing_length_a, 0.5, NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsGetPositiveScalar("-mixing_length_peak_amplitude", &P->mixing_length_b, 0.5, NULL);
    CHKERRQ(ierr);
  }
  else if (P->mixing_length == 2)
  {
    /* if P->mixing_length==2 then mixing length is constant, and these values are not required */
    P->mixing_length_a = 0.0;
    P->mixing_length_b = 0.0;
  }
  else
  {
    SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "-mixing_length must be 1 or 2 (not %d)", P->mixing_length);
  }

  /* option to include a mid-mantle layer
     this is non-dimensional fractional radius, as with coresize
     note this only influences the radial mixing length (not viscosity etc.) */
  {
    PetscBool set;
    ierr = PetscOptionsGetPositiveScalar("-layer_interface_radius", &P->layer_interface_radius, P->coresize, &set);
    CHKERRQ(ierr);
    if (set)
    {
      if ((P->layer_interface_radius < P->coresize || P->layer_interface_radius > 1.0))
      {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "-layer_interface_radius %f must be greater than -coresize and less than 1.0", P->layer_interface_radius);
      }
    }
  }

  P->Mg_Si0 = 0.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-Mg_Si0", &P->Mg_Si0, NULL);
  CHKERRQ(ierr);
  P->Mg_Si1 = 0.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-Mg_Si1", &P->Mg_Si1, NULL);
  CHKERRQ(ierr);

  /* initial condition for interior */
  ierr = PetscOptionsGetPositiveInt("-IC_INTERIOR", &P->IC_INTERIOR, 1, NULL);
  CHKERRQ(ierr);

  ierr = PetscStrcpy(P->ic_interior_filename, "restart.json");
  CHKERRQ(ierr);

  if (P->IC_INTERIOR == 2)
  {
    ierr = PetscOptionsGetString(NULL, NULL, "-ic_interior_filename", P->ic_interior_filename, PETSC_MAX_PATH_LEN, NULL);
    CHKERRQ(ierr);
  }

  /* start time (years) */
  P->t0 = 0.0;
  PetscScalar t0 = 0.0;
  PetscBool t0_set = PETSC_FALSE;
  ierr = PetscOptionsGetReal(NULL, NULL, "-t0", &t0, &t0_set);
  CHKERRQ(ierr);
  if (t0_set)
    P->t0 = t0;
  P->t0 /= SC->TIMEYRS; // non-dimensional for time stepping

  /* P->t0 could be subsequently overwritten by set_ic_interior() */

  P->ic_melt_pressure = 30.0; // GPa
  if (P->IC_INTERIOR == 3)
  {
    ierr = PetscOptionsGetScalar(NULL, NULL, "-ic_melt_pressure", &P->ic_melt_pressure, NULL);
    CHKERRQ(ierr);
  }

  /* initial temperature at top of adiabat (K) */
  P->ic_adiabat_temperature = 2000.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-ic_adiabat_temperature", &P->ic_adiabat_temperature, NULL);
  CHKERRQ(ierr);
  P->ic_adiabat_temperature /= SC->TEMP;

  /* initial temperature gradient (K/m) */
  P->ic_dTdr = -3.0e-04;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-ic_dTdr", &P->ic_dTdr, NULL);
  CHKERRQ(ierr);
  P->ic_dTdr /= SC->DTDR;

  P->ic_steady_state_energy = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-ic_steady_state_energy", &P->ic_steady_state_energy, NULL);
  CHKERRQ(ierr);

  /* initial temperature at the surface */
  P->ic_surface_temperature = -1;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-ic_surface_temperature", &P->ic_surface_temperature, NULL);
  CHKERRQ(ierr);
  if (P->ic_surface_temperature > 0.0)
  {
    P->ic_surface_temperature /= SC->TEMP;
  }

  /* initial temperature at the core-mantle boundary */
  P->ic_core_temperature = -1;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-ic_core_temperature", &P->ic_core_temperature, NULL);
  CHKERRQ(ierr);
  if (P->ic_core_temperature > 0.0)
  {
    P->ic_core_temperature /= SC->TEMP;
  }

  /* eos for determining mapping between radius and mass coordinate */
  ierr = EOSCreate(&P->eos_mesh, SPIDER_EOS_ADAMSWILLIAMSON);
  CHKERRQ(ierr);
  ierr = EOSSetUpFromOptions(P->eos_mesh, "adams_williamson", FC, SC);
  CHKERRQ(ierr);

  /* grain size (m) */
  P->grain = 1.0E-3;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-grain", &P->grain, NULL);
  CHKERRQ(ierr);
  P->grain /= SC->RADIUS;

  /* only allow energy transport upwards due to melt/solid separation if the cell below
     the cell interface contains melt/solid.  Otherwise, you can (unphysically) cool
     a cell due to melt migration even though it contains no melt.  For middle-out or
     more general scenarios, this requires more thought. */
  P->JGRAV_BOTTOM_UP = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-JGRAV_BOTTOM_UP", &P->JGRAV_BOTTOM_UP, NULL);
  CHKERRQ(ierr);

  {
    PetscInt CORE_BC = 0;
    PetscBool CORE_BCset = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL, NULL, "-CORE_BC", &CORE_BC, &CORE_BCset);
    CHKERRQ(ierr);
    if (CORE_BCset)
      P->CORE_BC = CORE_BC;
  }
  P->core_bc_value = 0.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-core_bc_value", &P->core_bc_value, NULL);
  CHKERRQ(ierr);

  /* gravity (m/s^2), must be negative */
  P->gravity = -10.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-gravity", &P->gravity, NULL);
  CHKERRQ(ierr);
  P->gravity /= SC->GRAVITY;

  /* melt fraction threshold for rheology */
  P->phi_critical = 0.4; // non dimensional
  ierr = PetscOptionsGetScalar(NULL, NULL, "-phi_critical", &P->phi_critical, NULL);
  CHKERRQ(ierr);
  /* melt fraction transition width for rheology */
  P->phi_width = 0.15; // non dimensional
  ierr = PetscOptionsGetScalar(NULL, NULL, "-phi_width", &P->phi_width, NULL);
  CHKERRQ(ierr);

  /* entropy of fusion */
  P->entropy_of_fusion = 900.0; // J/kg/K from eyeballing Bower et al. (2018), Fig. 3a.
  ierr = PetscOptionsGetScalar(NULL, NULL, "-entropy_of_fusion", &P->entropy_of_fusion, NULL);
  CHKERRQ(ierr);
  P->entropy_of_fusion /= SC->ENTROPY;

  /* core density (kg/m^3) */
  P->rho_core = 10738.332568062382;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-rho_core", &P->rho_core, NULL);
  CHKERRQ(ierr);
  P->rho_core /= SC->DENSITY;

  /* core mass (calculated, non-dimensional scaled mass) */
  P->coremass = 1.0 / 3.0 * PetscPowScalar(P->coresize, 3.0) * PetscPowScalar(P->radius, 3.0);
  P->coremass *= P->rho_core;

  /* heat capacity of core (J/kg/K) */
  P->cp_core = 880.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-cp_core", &P->cp_core, NULL);
  CHKERRQ(ierr);
  P->cp_core /= SC->ENTROPY;

  /* mass-weighted average core temperature as a fraction */
  /* of CMB temperature (non-dimensional) */
  P->tfac_core_avg = 1.147;

  /* viscosity cut-offs */
  P->log10visc_min = -1.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-log10visc_min", &P->log10visc_min, NULL);
  CHKERRQ(ierr);
  if (P->log10visc_min > 0.0)
  {
    P->log10visc_min -= SC->LOG10VISC;
  }
  P->log10visc_max = -1.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-log10visc_max", &P->log10visc_max, NULL);
  CHKERRQ(ierr);
  if (P->log10visc_max > 0.0)
  {
    P->log10visc_max -= SC->LOG10VISC;
  }

  /* option to scale eddy diffusivities for temperature and chemistry, or set as constants */
  P->eddy_diffusivity_thermal = 1.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-eddy_diffusivity_thermal", &P->eddy_diffusivity_thermal, NULL);
  CHKERRQ(ierr);
  /* if input is negative, then set as constant.  Retain negative sign to use as flag */
  if (P->eddy_diffusivity_thermal < 0.0)
  {
    /* must scale */
    P->eddy_diffusivity_thermal /= SC->KAPPA;
  }
  /* otherwise, we just scale the calculated eddy diffusivity by the user-specified constant */

  P->eddy_diffusivity_chemical = 1.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-eddy_diffusivity_chemical", &P->eddy_diffusivity_chemical, NULL);
  CHKERRQ(ierr);
  /* if input is negative, then set as constant.  Retain negative sign to use as flag */
  if (P->eddy_diffusivity_chemical < 0.0)
  {
    /* must scale */
    P->eddy_diffusivity_chemical /= SC->KAPPA;
  }
  /* otherwise, we just scale the calculated eddy diffusivity by the user-specified constant */

  /* viscous lid added by Rob Spaargaren */
  P->VISCOUS_LID = 0;
  ierr = PetscOptionsGetInt(NULL, NULL, "-VISCOUS_LID", &P->VISCOUS_LID, NULL);
  CHKERRQ(ierr);
  P->lid_log10visc = 0.0;
  P->lid_thickness = 0.0; /* metres */
  if (P->VISCOUS_LID)
  {
    ierr = PetscOptionsGetScalar(NULL, NULL, "-lid_log10visc", &P->lid_log10visc, NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(NULL, NULL, "-lid_thickness", &P->lid_thickness, NULL);
    CHKERRQ(ierr);
    P->lid_thickness /= SC->RADIUS;
  }

  /* core boundary condition */
  P->CORE_BC = MO_CORE_TYPE_COOLING;
  {
    PetscInt CORE_BC = 0;
    PetscBool CORE_BCset = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL, NULL, "-CORE_BC", &CORE_BC, &CORE_BCset);
    CHKERRQ(ierr);
    if (CORE_BCset)
      P->CORE_BC = CORE_BC;
  }
  P->core_bc_value = 0.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-core_bc_value", &P->core_bc_value, NULL);
  CHKERRQ(ierr);
  /* this switch is kept for legacy purposes, but the universal boundary condition for the CMB
     could in principle (in the future) be controlled by a single parameter:
         P->core_bc_value = -1 would be isothermal (Ecmb = Ecore)
         P->core_bc_value = 0 would be simple core cooling (no core heat flux)
         P->core_bc_value > 0 would be prescribed heat flux from core into mantle */
  switch (P->CORE_BC)
  {
  case 1:
    // CORE_BC = MO_CORE_TYPE_COOLING: simple core cooling
    P->core_bc_value = 0.0;
    break;
  case 2:
    // CORE_BC = MO_CORE_TYPE_HEAT_FLUX: heat flux (prescribed)
    P->core_bc_value /= SC->FLUX;
    break;
  case 3:
    /* CORE_BC = MO_CORE_TYPE_TEMPERATURE: temperature
       do nothing */
    break;
  default:
    SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unsupported CORE_BC value %d provided", P->CORE_BC);
    break;
  }

  /* Look for command-line option to determine number of radionuclides
     and options prefix for each e.g. -radionuclides_names al26, k40 */
  P->n_radionuclides = 0;
  {
    char *prefixes[SPIDER_MAX_RADIONUCLIDES];
    PetscInt n_radionuclides = SPIDER_MAX_RADIONUCLIDES;
    PetscBool set;

    ierr = PetscOptionsGetStringArray(NULL, NULL, "-radionuclide_names", prefixes, &n_radionuclides, &set);
    CHKERRQ(ierr);
    if (set)
    {
      PetscInt r;

      P->n_radionuclides = n_radionuclides;
      for (r = 0; r < P->n_radionuclides; ++r)
      {
        ierr = RadionuclideParametersCreate(&P->radionuclide_parameters[r]);
        CHKERRQ(ierr);
        ierr = PetscStrncpy(P->radionuclide_parameters[r]->prefix, prefixes[r], sizeof(P->radionuclide_parameters[r]->prefix));
        CHKERRQ(ierr);
        ierr = PetscFree(prefixes[r]);
        CHKERRQ(ierr);
        ierr = RadionuclideParametersSetFromOptions(P->radionuclide_parameters[r], SC);
        CHKERRQ(ierr);
      }
    }
  }

  /* Look for command-line option to determine number of phases
     and options prefix for each e.g. -phase_names melt,solid */
  P->n_phases = 0;
  {
    char *prefixes[SPIDER_MAX_PHASES];
    PetscInt n_phases = SPIDER_MAX_PHASES;
    PetscBool set;

    ierr = PetscOptionsGetStringArray(NULL, NULL, "-phase_names", prefixes, &n_phases, &set);
    CHKERRQ(ierr);
    /* If there is a single phase, use this "pure" EOS. If there are two phases,
       form a composite, assuming melt,solid ordering */
    if (set)
    {
      P->n_phases = n_phases;
      if (P->n_phases > SPIDER_MAX_PHASES)
        SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP, "%D phases specified, but the maximum is %D", P->n_phases, SPIDER_MAX_PHASES);
      for (PetscInt r = 0; r < P->n_phases; ++r)
      {
        char buf[1024]; /* max size */
        PetscInt type = 1;

        ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefixes[r], "_TYPE");
        CHKERRQ(ierr);
        ierr = PetscOptionsGetInt(NULL, NULL, buf, &type, NULL);
        CHKERRQ(ierr);
        switch (type)
        {
        case 1:
          ierr = EOSCreate(&P->eos_phases[r], SPIDER_EOS_LOOKUP);
          CHKERRQ(ierr);
          break;
        default:
          SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unrecognized type code %D", type);
        }
        ierr = EOSSetUpFromOptions(P->eos_phases[r], prefixes[r], FC, SC);
        CHKERRQ(ierr);
        ierr = PetscFree(prefixes[r]);
        CHKERRQ(ierr);
      }
      if (P->n_phases == 1)
      {
        P->eos = P->eos_phases[0];
      }
      else if (P->n_phases == 2)
      {
        ierr = EOSCreate(&P->eos, SPIDER_EOS_COMPOSITE);
        CHKERRQ(ierr);
        ierr = EOSCompositeSetSubEOS(P->eos, P->eos_phases, P->n_phases);
        CHKERRQ(ierr);
        ierr = EOSSetUpFromOptions(P->eos, "composite", FC, SC);
        CHKERRQ(ierr);
      }
      else
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Only one or two phases are supported");
    }
    else
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "You must supply the -phase_names option");
  }

  /* Energy terms to include */
  P->CONDUCTION = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-CONDUCTION", &P->CONDUCTION, NULL);
  CHKERRQ(ierr);
  P->CONVECTION = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-CONVECTION", &P->CONVECTION, NULL);
  CHKERRQ(ierr);
  P->HTIDAL = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-HTIDAL", &P->HTIDAL, NULL);
  CHKERRQ(ierr);
  P->MIXING = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-MIXING", &P->MIXING, NULL);
  CHKERRQ(ierr);
  P->SEPARATION = 1;
  ierr = PetscOptionsGetInt(NULL, NULL, "-SEPARATION", &P->SEPARATION, NULL);
  CHKERRQ(ierr);

  /* separation and mixing only relevant for multiphase systems */
  if (P->n_phases == 1)
  {
    P->MIXING = PETSC_FALSE;
    P->SEPARATION = 0;
  }

  ierr = AtmosphereParametersSetFromOptions(P, SC, FC);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode AtmosphereParametersSetFromOptions(Parameters P, const ScalingConstants SC, const FundamentalConstants FC)
{
  PetscErrorCode ierr;
  AtmosphereParameters Ap = P->atmosphere_parameters;
  PetscInt v;

  PetscFunctionBeginUser;

  /* initial condition for atmosphere */
  Ap->IC_ATMOSPHERE = 1;
  ierr = PetscOptionsGetInt(NULL, NULL, "-IC_ATMOSPHERE", &Ap->IC_ATMOSPHERE, NULL);
  CHKERRQ(ierr);

  ierr = PetscStrcpy(Ap->ic_atmosphere_filename, "restart.json");
  CHKERRQ(ierr);
  if ((Ap->IC_ATMOSPHERE == 2) || (Ap->IC_ATMOSPHERE == 3))
  {
    ierr = PetscOptionsGetString(NULL, NULL, "-ic_atmosphere_filename", Ap->ic_atmosphere_filename, PETSC_MAX_PATH_LEN, NULL);
    CHKERRQ(ierr);
  }

  /* solve for surface radiative balance during time-stepping */
  Ap->SURFACE_BC_ACC = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-SURFACE_BC_ACC", &Ap->SURFACE_BC_ACC, NULL);
  CHKERRQ(ierr);

  /* (top) surface boundary condition */
  Ap->SURFACE_BC = MO_ATMOSPHERE_TYPE_GREY_BODY;
  {
    PetscInt SURFACE_BC = 0;
    PetscBool SURFACE_BCset = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL, NULL, "-SURFACE_BC", &SURFACE_BC, &SURFACE_BCset);
    CHKERRQ(ierr);
    if (SURFACE_BCset)
      Ap->SURFACE_BC = SURFACE_BC;
  }
  Ap->surface_bc_value = 0.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-surface_bc_value", &Ap->surface_bc_value, NULL);
  CHKERRQ(ierr);
  switch (Ap->SURFACE_BC)
  {
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
       using plane-parallel radiative equilibrium model of Abe and Matsui (1985)
       do nothing */
    break;
  case 4:
    /* MO_ATMOSPHERE_TYPE_HEAT_FLUX: heat flux (prescribed) */
    Ap->surface_bc_value /= SC->FLUX;
    break;
  case 5:
    /* SURFACE_BC = MO_ATMOSPHERE_TYPE_ENTROPY: entropy
       do nothing */
    break;
  default:
    SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unsupported SURFACE_BC value %d provided", Ap->SURFACE_BC);
    break;
  }

  /* use psuedo volatiles */
  Ap->PSEUDO_VOLATILES = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-PSEUDO_VOLATILES", &Ap->PSEUDO_VOLATILES, NULL);
  CHKERRQ(ierr);

  /* emissivity is constant for SURFACE_BC != MO_ATMOSPHERE_TYPE_VOLATILES */
  Ap->emissivity0 = 1.0; // non-dimensional
  ierr = PetscOptionsGetScalar(NULL, NULL, "-emissivity0", &Ap->emissivity0, NULL);
  CHKERRQ(ierr);

  Ap->THERMAL_ESCAPE = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-THERMAL_ESCAPE", &Ap->THERMAL_ESCAPE, NULL);
  CHKERRQ(ierr);

  Ap->CONSTANT_ESCAPE = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-CONSTANT_ESCAPE", &Ap->CONSTANT_ESCAPE, NULL);
  CHKERRQ(ierr);

  /* Zahnle et al. (2019), Eqn 3, accounting for both diffusion limited and
     energy limited flux */
  Ap->ZAHNLE_ESCAPE = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-ZAHNLE_ESCAPE", &Ap->ZAHNLE_ESCAPE, NULL);
  CHKERRQ(ierr);

  /* for Zahnle et al. (2019), Eqn. 3, need S(t) = F_xuv / F_xuv_present (see Eqn 2) */
  Ap->solar_xuv_factor = 1.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-solar_xuv_factor", &Ap->solar_xuv_factor, NULL);
  CHKERRQ(ierr);

  /* equilibrium temperature of the planet (K) */
  Ap->teqm = 273.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-teqm", &Ap->teqm, NULL);
  CHKERRQ(ierr);
  Ap->teqm /= SC->TEMP;

  Ap->tsurf_poststep_change = -1; // (K) (negative value is OFF)
  ierr = PetscOptionsGetScalar(NULL, NULL, "-tsurf_poststep_change", &Ap->tsurf_poststep_change, NULL);
  CHKERRQ(ierr);
  Ap->tsurf_poststep_change /= SC->TEMP;

  /* for radiative boundary condition at the top surface
     dT = param_utbl_const * [Surface temperature]**3 */
  Ap->PARAM_UTBL = PETSC_TRUE;
  Ap->param_utbl_const = 1.0e-7;
  ierr = PetscOptionsGetBool(NULL, NULL, "-PARAM_UTBL", &Ap->PARAM_UTBL, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL, NULL, "-param_utbl_const", &Ap->param_utbl_const, NULL);
  CHKERRQ(ierr);
  if (Ap->PARAM_UTBL)
  {
    Ap->param_utbl_const *= PetscSqr(SC->TEMP);
  }
  else
  {
    Ap->param_utbl_const = 0.0;
  }

  /* below here are only used for SURFACE_BC = MO_ATMOSPHERE_TYPE_VOLATILES */

  /* atmosphere reference pressure (Pa) */
  Ap->P0 = 101325.0; // Pa (= 1 atm)
  ierr = PetscOptionsGetScalar(NULL, NULL, "-P0", &Ap->P0, NULL);
  CHKERRQ(ierr);
  Ap->P0 /= SC->PRESSURE;

  Ap->gravity_ptr = &P->gravity;
  Ap->radius_ptr = &P->radius;
  Ap->VOLATILE_ptr = &SC->VOLATILE;

  /* Look for command-line option to determine number of volatiles
     and options prefix for each, e.g -volatile_names CO2,H2O */
  Ap->n_volatiles = 0;
  {
    char *prefixes[SPIDER_MAX_VOLATILE_SPECIES];
    PetscInt n_volatiles = SPIDER_MAX_VOLATILE_SPECIES;
    PetscBool set;

    ierr = PetscOptionsGetStringArray(NULL, NULL, "-volatile_names", prefixes, &n_volatiles, &set);
    CHKERRQ(ierr);
    if (set)
    {
      PetscInt v;
      Ap->n_volatiles = n_volatiles;
      for (v = 0; v < Ap->n_volatiles; ++v)
      {
        ierr = VolatileParametersCreate(&Ap->volatile_parameters[v]);
        CHKERRQ(ierr);
        ierr = PetscStrncpy(Ap->volatile_parameters[v]->prefix, prefixes[v], sizeof(Ap->volatile_parameters[v]->prefix));
        CHKERRQ(ierr);
        ierr = PetscFree(prefixes[v]);
        CHKERRQ(ierr);
      }
    }
  }

  /* Get command-line values for all volatiles species */
  for (v = 0; v < Ap->n_volatiles; ++v)
  {
    ierr = VolatileParametersSetFromOptions(Ap->volatile_parameters[v], Ap, SC, FC);
    CHKERRQ(ierr);
  }

  /* Reactions: look for command-line options to determine the number of reactions
     and options for each. These include named reactions and "simple" reactions.
     These are defined with respect to the prefixes for the volatiles, which are translated to integer ids */
  Ap->n_reactions = 0;

  /* Special "named" reactions */
  {
    PetscBool flg;

    ierr = PetscOptionsGetBool(NULL, NULL, "-reaction_ammonia1", NULL, &flg);
    CHKERRQ(ierr);
    if (flg)
    {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS)
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source", SPIDER_MAX_REACTIONS);
      ierr = ReactionParametersCreateAmmonia1(&Ap->reaction_parameters[Ap->n_reactions], Ap, SC);
      CHKERRQ(ierr);
      ++Ap->n_reactions;
    }

    /* ideally, should not allow the user to select both water reactions at the same time */
    ierr = PetscOptionsGetBool(NULL, NULL, "-reaction_carbondioxide_IVTANTHERMO", NULL, &flg);
    CHKERRQ(ierr);
    if (flg)
    {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS)
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source", SPIDER_MAX_REACTIONS);
      ierr = ReactionParametersCreateCarbonDioxideIVTANTHERMO(&Ap->reaction_parameters[Ap->n_reactions], Ap, SC);
      CHKERRQ(ierr);
      ++Ap->n_reactions;
    }

    /* ideally, should not allow the user to select both carbon dioxide reactions at the same time */
    ierr = PetscOptionsGetBool(NULL, NULL, "-reaction_carbondioxide_JANAF", NULL, &flg);
    CHKERRQ(ierr);
    if (flg)
    {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS)
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source", SPIDER_MAX_REACTIONS);
      ierr = ReactionParametersCreateCarbonDioxideJANAF(&Ap->reaction_parameters[Ap->n_reactions], Ap, SC);
      CHKERRQ(ierr);
      ++Ap->n_reactions;
    }

    ierr = PetscOptionsGetBool(NULL, NULL, "-reaction_methane_IVTANTHERMO", NULL, &flg);
    CHKERRQ(ierr);
    if (flg)
    {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS)
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source", SPIDER_MAX_REACTIONS);
      ierr = ReactionParametersCreateMethaneIVTANTHERMO(&Ap->reaction_parameters[Ap->n_reactions], Ap, SC);
      CHKERRQ(ierr);
      ++Ap->n_reactions;
    }

    /* ideally, should not allow the user to select both water reactions at the same time */
    ierr = PetscOptionsGetBool(NULL, NULL, "-reaction_water_IVTANTHERMO", NULL, &flg);
    CHKERRQ(ierr);
    if (flg)
    {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS)
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source", SPIDER_MAX_REACTIONS);
      ierr = ReactionParametersCreateWaterIVTANTHERMO(&Ap->reaction_parameters[Ap->n_reactions], Ap, SC);
      CHKERRQ(ierr);
      ++Ap->n_reactions;
    }

    ierr = PetscOptionsGetBool(NULL, NULL, "-reaction_water_JANAF", NULL, &flg);
    CHKERRQ(ierr);
    if (flg)
    {
      if (Ap->n_reactions >= SPIDER_MAX_REACTIONS)
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Too many reactions. Increase SPIDER_MAX_REACTIONS (currently %d) in the source", SPIDER_MAX_REACTIONS);
      ierr = ReactionParametersCreateWaterJANAF(&Ap->reaction_parameters[Ap->n_reactions], Ap, SC);
      CHKERRQ(ierr);
      ++Ap->n_reactions;
    }
  }

  Ap->OXYGEN_FUGACITY = OXYGEN_FUGACITY_NONE;
  {
    PetscInt OXYGEN_FUGACITY = 0;
    PetscBool OXYGEN_FUGACITYset = PETSC_FALSE;
    ierr = PetscOptionsGetInt(NULL, NULL, "-OXYGEN_FUGACITY", &OXYGEN_FUGACITY, &OXYGEN_FUGACITYset);
    CHKERRQ(ierr);
    /* only set OXYGEN_FUGACITY if there are reactions, since fO2 is
       used to also adjust the total atmosphere pressure */
    if (OXYGEN_FUGACITYset && Ap->n_reactions)
      Ap->OXYGEN_FUGACITY = OXYGEN_FUGACITY;
  }

  /* default offset, assuming IW buffer is 0.5 from Sossi et al. (2020) */
  Ap->OXYGEN_FUGACITY_offset = 0.5;
  ierr = PetscOptionsGetScalar(NULL, NULL, "-OXYGEN_FUGACITY_offset", &Ap->OXYGEN_FUGACITY_offset, NULL);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PrintParameters(Parameters const P)
{
  PetscErrorCode ierr;
  PetscInt i;
  ScalingConstants const SC = P->scaling_constants;
  AtmosphereParameters const Ap = P->atmosphere_parameters;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "**************** Magma Ocean | Parameters **************\n\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15s %s\n", "[Scaling]", "", "Value", "Units");
  CHKERRQ(ierr);
  /* these are primary scalings (can be user specified) */
  ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15.6g %-6s\n", "Radius", "", (double)SC->RADIUS, "m");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15.6g %-6s\n", "Temperature", "", (double)SC->TEMP, "K");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15.6g %-6s\n", "Entropy", "", (double)SC->ENTROPY, "J/kg/K");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15.6g %-6s\n", "Pressure", "", (double)SC->PRESSURE, "Pa");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15.6g %-6s\n", "Mass", "", (double)SC->MASS, "kg");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15.6g %-6s\n", "Volatile", "", (double)SC->VOLATILE, "mass fraction");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15.6g %-6s (%.6g years)\n", "Time", "", (double)SC->TIME, "s", (double)SC->TIMEYRS);
  CHKERRQ(ierr);
  /* next are derived from primary scalings and are useful for analysing the
     scaling of the numerical system of equations */
  // ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------\n"                                                 );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-28s %-2s %-15.6g %-6s\n", "dT/dr (temperature gradient)", "", (double)SC->DTDR, "K/m");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-28s %-2s %-15.6g %-6s\n", "T     (surface temperature)", "", (double)SC->TEMP, "K");
  CHKERRQ(ierr);
  PetscScalar const PRESSUREBAR = SC->PRESSURE * 1.0E-5; // bar usual units for atmosphere
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-28s %-2s %-15.6g %-6s (%.6g bar)\n", "Pv    (atmospheric pressure)", "", (double)SC->PRESSURE, "Pa", (double)PRESSUREBAR);
  CHKERRQ(ierr);
  PetscScalar const VOLMASS = SC->VOLATILE * SC->MASS; // volatile reservoir mass scaling
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-28s %-2s %-15.6g %-6s\n", "Mv    (volatile mass)", "", (double)VOLMASS, "kg");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15s %-15s %s\n", "[Parameter]", "Non-dim. Value", "Dim. Value", "Units");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15.6g %-15.6g %s\n", "dtmacro", (double)P->dtmacro, (double)(P->dtmacro * SC->TIME), "s");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15d\n", "nstepsmacro", P->nstepsmacro);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15d\n", "numpts_b", P->numpts_b);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15d\n", "numpts_s", P->numpts_s);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-15s %-15.6g %-15.6g %s\n", "ic_adiabat_temperature", (double)P->ic_adiabat_temperature, (double)(P->ic_adiabat_temperature * SC->TEMP), "K");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
  CHKERRQ(ierr);
  if (Ap->n_volatiles > 0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[Volatile] prefix/name\n");
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
    CHKERRQ(ierr);
    for (i = 0; i < Ap->n_volatiles; ++i)
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "%-10D %-15s (.. additional parameters omitted ..)\n", i, Ap->volatile_parameters[i]->prefix);
      CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
    CHKERRQ(ierr);
  }
  if (P->postStepActive)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "PostStep logic active\n");
    CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%-30s %s\n", "Output Directory", P->outputDirectory);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n********************************************************\n");
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
 ******************************************************************************
 * Create and Destroy functions for parameter structs
 ******************************************************************************
 */

/* Note: functions to Create and Destroy eos-related structs are in eos.c */
/* Note: functions to Create and Destroy ScalingConstants and FundamentalConstants are in constants.c */

static PetscErrorCode VolatileParametersCreate(VolatileParameters *volatile_parameters_ptr)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscMalloc1(1, volatile_parameters_ptr);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode VolatileParametersDestroy(VolatileParameters *volatile_parameters_ptr)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscFree(*volatile_parameters_ptr);
  CHKERRQ(ierr);
  *volatile_parameters_ptr = NULL;

  PetscFunctionReturn(0);
}

static PetscErrorCode RadionuclideParametersCreate(RadionuclideParameters *radionuclide_parameters_ptr)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscMalloc1(1, radionuclide_parameters_ptr);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode RadionuclideParametersDestroy(RadionuclideParameters *radionuclide_parameters_ptr)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscFree(*radionuclide_parameters_ptr);
  CHKERRQ(ierr);
  *radionuclide_parameters_ptr = NULL;

  PetscFunctionReturn(0);
}

static PetscErrorCode AtmosphereParametersCreate(AtmosphereParameters *atmosphere_parameters_ptr)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscMalloc1(1, atmosphere_parameters_ptr);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode AtmosphereParametersDestroy(AtmosphereParameters *atmosphere_parameters_ptr)
{
  PetscErrorCode ierr;
  PetscInt i;
  AtmosphereParameters Ap = *atmosphere_parameters_ptr;

  PetscFunctionBeginUser;

  for (i = 0; i < Ap->n_reactions; ++i)
  {
    ierr = ReactionParametersDestroy(&Ap->reaction_parameters[i]);
    CHKERRQ(ierr);
  }
  Ap->n_reactions = 0;

  for (i = 0; i < Ap->n_volatiles; ++i)
  {
    /* if pseudo-volatile, destroy Interp1d */
    if (Ap->PSEUDO_VOLATILES)
    {
      ierr = Interp1dDestroy(&Ap->volatile_parameters[i]->TP_interp);
      CHKERRQ(ierr);
    }
    ierr = VolatileParametersDestroy(&Ap->volatile_parameters[i]);
    CHKERRQ(ierr);
  }
  Ap->n_volatiles = 0;

  ierr = PetscFree(*atmosphere_parameters_ptr);
  CHKERRQ(ierr);
  *atmosphere_parameters_ptr = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode ParametersCreate(Parameters *parameters_ptr)
{
  PetscErrorCode ierr;
  Parameters P;

  /* main parameters struct */
  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, parameters_ptr);
  CHKERRQ(ierr);
  P = *parameters_ptr;

  /* nested structs */
  ierr = ScalingConstantsCreate(&P->scaling_constants);
  CHKERRQ(ierr);
  ierr = FundamentalConstantsCreate(&P->fundamental_constants);
  CHKERRQ(ierr);
  ierr = AtmosphereParametersCreate(&P->atmosphere_parameters);
  CHKERRQ(ierr);

  /* populate structs with data */
  ierr = ParametersSetFromOptions(P);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ParametersDestroy(Parameters *parameters_ptr)
{
  PetscErrorCode ierr;
  PetscInt i;
  Parameters P = *parameters_ptr;

  PetscFunctionBeginUser;

  ierr = ScalingConstantsDestroy(&P->scaling_constants);
  CHKERRQ(ierr);
  ierr = FundamentalConstantsDestroy(&P->fundamental_constants);
  CHKERRQ(ierr);
  ierr = AtmosphereParametersDestroy(&P->atmosphere_parameters);
  CHKERRQ(ierr);

  /* radionuclides */
  for (i = 0; i < P->n_radionuclides; ++i)
  {
    ierr = RadionuclideParametersDestroy(&P->radionuclide_parameters[i]);
    CHKERRQ(ierr);
  }
  P->n_radionuclides = 0;

  /* radius to mass coordinate eos */
  ierr = EOSDestroy(&P->eos_mesh);

  /* EOS / phases */
  for (i = 0; i < P->n_phases; ++i)
  {
    ierr = EOSDestroy(&P->eos_phases[i]);
    CHKERRQ(ierr);
  }
  /* destroy composite */
  if (P->n_phases == 2)
  {
    ierr = EOSDestroy(&P->eos);
    CHKERRQ(ierr);
  }

  P->n_phases = 0;
  P->eos = NULL;

  ierr = PetscFree(*parameters_ptr);
  CHKERRQ(ierr);
  *parameters_ptr = NULL;
  PetscFunctionReturn(0);
}
