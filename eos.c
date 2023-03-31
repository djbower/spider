#include "eos.h"
#include "util.h"

/* Prototypes for helper functions used in interface functions */
static PetscScalar GetCompositionalViscosityPrefactor(PetscScalar);
static PetscErrorCode LookupFilenameSet(const char *, const char *, char *, PetscBool *);

/* EOS interface functions (public API) */
PetscErrorCode EOSCheckType(EOS eos, EOSType type, PetscBool *same)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscStrcmp(type, eos->type, same);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EOSCreate(EOS *p_eos, EOSType type)
{
  PetscErrorCode ierr;
  PetscBool flg;
  EOS eos;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, p_eos);
  CHKERRQ(ierr);
  eos = *p_eos;
  eos->type = type;
  eos->is_setup = PETSC_FALSE;
  ierr = PetscStrcmp(type, SPIDER_EOS_LOOKUP, &flg);
  CHKERRQ(ierr);
  if (flg)
  {
    ierr = EOSCreate_Lookup(eos);
    CHKERRQ(ierr);
  }
  ierr = PetscStrcmp(type, SPIDER_EOS_COMPOSITE, &flg);
  CHKERRQ(ierr);
  if (flg)
  {
    ierr = EOSCreate_Composite(eos);
    CHKERRQ(ierr);
  }
  ierr = PetscStrcmp(type, SPIDER_EOS_ADAMSWILLIAMSON, &flg);
  CHKERRQ(ierr);
  if (flg)
  {
    ierr = EOSCreate_AdamsWilliamson(eos);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode EOSEval(const EOS eos, PetscScalar P, PetscScalar S, EOSEvalData *eval)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  if (!eos->is_setup)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "You must cannot evaluate the EOS before setting it up");

  /* Implementation-specific logic */
  if (!eos->eval)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Incomplete EOS object: no implementation for EOSEval has been provided");
  ierr = (*(eos->eval))(eos, P, S, eval);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode EOSDestroy(EOS *p_eos)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  if (p_eos)
  {
    EOS eos = *p_eos;

    ierr = (*eos->destroy)(eos);
    CHKERRQ(ierr);
    if (eos->PHASE_BOUNDARY)
    {
      ierr = Interp1dDestroy(&eos->phase_boundary);
      CHKERRQ(ierr);
    }
    ierr = PetscFree(*p_eos);
    CHKERRQ(ierr);
    *p_eos = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EOSSetUpFromOptions(EOS eos, const char *prefix, const FundamentalConstants FC, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  char buf[1024]; /* max size */

  PetscFunctionBegin;
  if (eos->setupfromoptions)
  {
    ierr = (*eos->setupfromoptions)(eos, prefix, FC, SC);
    CHKERRQ(ierr);
  }

  /* conductivity (w/m/K) */
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_cond");
  CHKERRQ(ierr);
  eos->cond = 4.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->cond, NULL);
  CHKERRQ(ierr);
  eos->cond /= SC->COND;

  /* viscosity-related, may eventually move into their own struct */
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_log10visc");
  CHKERRQ(ierr);
  eos->log10visc = 21.0; // note default is for solid only
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->log10visc, NULL);
  CHKERRQ(ierr);
  eos->log10visc -= SC->LOG10VISC;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_activation_energy");
  CHKERRQ(ierr);
  eos->activation_energy = 0.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->activation_energy, NULL);
  CHKERRQ(ierr);
  /* scalings for Ea and Va include the gas constant that appears in the
     denominator of the Arrhenius viscosity law, so we then do not have
     to pass FC->GAS to the viscosity functions */
  eos->activation_energy /= SC->ENERGY * FC->GAS;

  /* activation volume (m^3/mol) */
  /* The numerical value in units of m^3/mol is the same as that in units of J/mol/Pa */
  /* You can convince yourself of this by using the scalings for ENERGY and PRESSURE to
     see that this is true */
  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_activation_volume");
  CHKERRQ(ierr);
  eos->activation_volume = 0.0;
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->activation_volume, NULL);
  CHKERRQ(ierr);
  /* as with activation energy, include gas constant in denominator */
  eos->activation_volume *= SC->PRESSURE / (SC->ENERGY * FC->GAS);

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_activation_volume_pressure_scale");
  CHKERRQ(ierr);
  eos->activation_volume_pressure_scale = -1.0; /* negative is not set */
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->activation_volume_pressure_scale, NULL);
  CHKERRQ(ierr);
  eos->activation_volume_pressure_scale /= SC->PRESSURE;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_visc_comp");
  CHKERRQ(ierr);
  eos->visc_comp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->visc_comp, NULL);
  CHKERRQ(ierr);
  /* no scaling necessary since this is a ratio */

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_visc_ref_temp");
  CHKERRQ(ierr);
  eos->visc_ref_temp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->visc_ref_temp, NULL);
  CHKERRQ(ierr);
  eos->visc_ref_temp /= SC->TEMP;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_visc_ref_pressure");
  CHKERRQ(ierr);
  eos->visc_ref_pressure = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->visc_ref_pressure, NULL);
  CHKERRQ(ierr);
  eos->visc_ref_pressure /= SC->PRESSURE;

  ierr = PetscSNPrintf(buf, sizeof(buf), "%s%s%s", "-", prefix, "_visc_ref_comp");
  CHKERRQ(ierr);
  eos->visc_ref_comp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL, NULL, buf, &eos->visc_ref_comp, NULL);
  CHKERRQ(ierr);
  /* no scaling necessary since this is a ratio */

  /* phase boundary */
  ierr = LookupFilenameSet("_phase_boundary", prefix, eos->phase_boundary_filename, &eos->PHASE_BOUNDARY);
  CHKERRQ(ierr);
  if (eos->PHASE_BOUNDARY)
  {
    ierr = Interp1dCreateAndSet(eos->phase_boundary_filename, &eos->phase_boundary, SC->PRESSURE, SC->ENTROPY);
    CHKERRQ(ierr);
  }

  /* Mark as set up */
  eos->is_setup = PETSC_TRUE;

  PetscFunctionReturn(0);
}

PetscErrorCode EOSGetPhaseBoundary(EOS eos, PetscScalar P, PetscScalar *boundary, PetscScalar *dboundary)
{
  PetscFunctionBeginUser;
  SetInterp1dValue(eos->phase_boundary, P, boundary, dboundary); /* entropy S and derivative dS/dP */
  PetscFunctionReturn(0);
}

PetscErrorCode EOSGetType(EOS eos, EOSType *type)
{
  PetscFunctionBeginUser;
  *type = eos->type;
  PetscFunctionReturn(0);
}

PetscErrorCode EOSEvalSetViscosity(EOS eos, EOSEvalData *eval)
{
  PetscScalar A, log10C, dP, dT;
  PetscScalar fac1 = 1.0, fac2 = 1.0;

  PetscFunctionBeginUser;

  /* reference viscosity */
  eval->log10visc = eos->log10visc; // i.e., log10(eta_0)

  /* temperature and pressure contribution
     A(T,P) = (E_a + V_a P) / RT - (E_a + V_a Pref)/ RTref
     eta = eta_0 * exp(A)
     log10(eta) = log10(eta0) + log10(exp(A))
     log10(eta) = P->eos2_parameters.log10visc + A/ln(10) */

  /* with Ps = activation_volume_pressure_scale:
         V_a(P) = V_a exp (-P/Ps) */

  A = 0.0;

  /* pin viscosity profile to reference values */
  if ((eos->visc_ref_pressure >= 0.0) || (eos->visc_ref_temp >= 0.0))
  {
    dT = (eos->visc_ref_temp - eval->T) / eos->visc_ref_temp;
    if (eos->activation_volume_pressure_scale > 0.0)
    {
      fac1 = PetscExpReal(-eval->P / eos->activation_volume_pressure_scale);
      fac2 = PetscExpReal(-eos->visc_ref_pressure / eos->activation_volume_pressure_scale);
    }
    /* else fac1 and fac2 retain unity scalings according to initialisation above */
    dP = fac1 * eval->P - fac2 * eos->visc_ref_pressure * (eval->T / eos->visc_ref_temp);
  }
  /* do not pin viscosity profile */
  else
  {
    dT = 1.0;
    dP = eval->P;
    if (eos->activation_volume_pressure_scale > 0.0)
    {
      dP *= PetscExpReal(-eval->P / eos->activation_volume_pressure_scale);
    }
  }

  if (eos->activation_energy > 0.0)
  {
    A += eos->activation_energy * dT;
  }
  if (eos->activation_volume > 0.0)
  {
    A += eos->activation_volume * dP;
  }

  /* division by R (gas constant) was already done during the scaling of parameters */
  A *= 1.0 / eval->T;
  eval->log10visc += A / PetscLogReal(10.0);

  /* compositional (Mg/Si) contribution */
  /* always pinned to some reference given by eos->visc_ref_comp */
  if (eos->visc_comp > 0.0)
  {
    log10C = GetCompositionalViscosityPrefactor(eos->visc_comp);
    log10C -= GetCompositionalViscosityPrefactor(eos->visc_ref_comp);
    eval->log10visc += log10C;
  }

  /* add viscous lid */

  /* add viscosity cutoff */

  PetscFunctionReturn(0);
}

/* Helper Functions */
static PetscScalar GetCompositionalViscosityPrefactor(PetscScalar Mg_Si)
{

  /* These expressions were worked out by Rob Spaargaren as part
     of his MSc thesis (2018) are are explained in Spaargaren et al. (2020) */

  /* Mg_Si is molar mantle Mg/Si */

  PetscScalar fac;

  if (Mg_Si <= 1.0)
    fac = 0.5185 * (1 - Mg_Si) / 0.3; //
  else if (Mg_Si <= 1.25)
    /* Earth has Mg/Si = 1.08 */
    fac = -1.4815 * (Mg_Si - 1) / 0.25; // -1.4815 = log10(0.033)
  else if (Mg_Si <= 1.5)
    fac = -2 + (0.5185) * (1.5 - Mg_Si) / 0.25; // 0.5185 = log10(0.033) - -2
  else
    /* Fp-rich composition (Ballmer et al. 2017) */
    fac = -2;

  return fac;
}

static PetscErrorCode LookupFilenameSet(const char *property, const char *prefix, char *lookup_filename, PetscBool *IS_SET)
{
  PetscErrorCode ierr;
  char buf1[1024]; /* max size */
  char buf2[1024]; /* max size */
  PetscBool set_rel_to_src, set;

  PetscFunctionBeginUser;

  /* Based on input options, determine which files to load.  Options ending
     with _rel_to_src indicate a path relative to the source code. In this
     case we prepend a string, SPIDER_ROOT_DIR_STR, and /. The corresponding
     option without this overrides. */

  /* check for relative path name */
  ierr = PetscSNPrintf(buf1, sizeof(buf1), "%s%s%s%s", "-", prefix, property, "_filename_rel_to_src");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, buf1, lookup_filename, PETSC_MAX_PATH_LEN, &set_rel_to_src);
  CHKERRQ(ierr);
  if (set_rel_to_src)
  {
    ierr = MakeRelativeToSourcePathAbsolute(lookup_filename);
    CHKERRQ(ierr);
  }
  /* check for absolute path name */
  ierr = PetscSNPrintf(buf2, sizeof(buf2), "%s%s%s%s", "-", prefix, property, "_filename");
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, buf2, lookup_filename, PETSC_MAX_PATH_LEN, &set);
  CHKERRQ(ierr);

  /* if IS_SET is NULL, then we require a valid lookup_filename to be returned */
  if (IS_SET == NULL)
  {
    /* must return a valid lookup_filename */
    if (!set && !set_rel_to_src)
    {
      SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "Missing argument %s or %s", buf1, buf2);
    }
  }

  /* if IS_SET is not NULL, then a valid lookup_filename is optional */
  if (IS_SET != NULL)
  {
    if (set || set_rel_to_src)
    {
      *IS_SET = PETSC_TRUE;
    }
    else
    {
      *IS_SET = PETSC_FALSE;
    }
  }

  /* absolute path name always overrides relative path name */
  if (set && set_rel_to_src)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%s%s%s%s%s", "Warning: ", buf1, " ignored because ", buf2, " provided\n");
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
