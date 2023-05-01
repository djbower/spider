#include "reaction.h"

static PetscScalar get_psurf_exponent(const ReactionParameters);
static PetscScalar get_reaction_quotient(const ReactionParameters, const Atmosphere *, PetscInt);
static PetscScalar get_reaction_quotient_time_derivative(const ReactionParameters, const Atmosphere *, const AtmosphereParameters, PetscInt);
static PetscScalar get_log10_equilibrium_constant(const ReactionParameters, PetscScalar);
static PetscScalar get_dlog10KdT(const ReactionParameters, PetscScalar);

/* Note: this could logically be included in parameters.c, but that file was getting crowded */

/* A named reaction */
PetscErrorCode ReactionParametersCreateMethaneIVTANTHERMO(ReactionParameters *reaction_parameters_ptr, const AtmosphereParameters Ap, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  PetscInt i, v;
  PetscBool flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, reaction_parameters_ptr);
  CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "methaneIVTANTHERMO";
  reaction_parameters->n_volatiles = 3;
  ierr = PetscMalloc3(2, &reaction_parameters->Keq_coeffs, reaction_parameters->n_volatiles, &reaction_parameters->stoichiometry, reaction_parameters->n_volatiles, &reaction_parameters->volatiles);
  CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0; // CO2
  reaction_parameters->stoichiometry[1] = -2.0; // H2
  reaction_parameters->stoichiometry[2] = 1.0;  // CH4
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = -16276;    // K
  reaction_parameters->Keq_coeffs[0] /= SC->TEMP; // non-dimensional
  reaction_parameters->Keq_coeffs[1] = -5.4738;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 1.0;

  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    reaction_parameters->volatiles[i] = -1; /* error value */
  for (v = 0; v < Ap->n_volatiles; ++v)
  {
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "CO2", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "H2", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[1] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "CH4", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[2] = v;
  }
  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    if (reaction_parameters->volatiles[i] == -1)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Didn't find required volatiles for reaction %s", reaction_parameters->type);
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersCreateAmmonia1(ReactionParameters *reaction_parameters_ptr, const AtmosphereParameters Ap, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  PetscInt i, v;
  PetscBool flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, reaction_parameters_ptr);
  CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "ammonia1";
  reaction_parameters->n_volatiles = 3;
  ierr = PetscMalloc3(2, &reaction_parameters->Keq_coeffs, reaction_parameters->n_volatiles, &reaction_parameters->stoichiometry, reaction_parameters->n_volatiles, &reaction_parameters->volatiles);
  CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0; // N2
  reaction_parameters->stoichiometry[1] = -3.0; // H2
  reaction_parameters->stoichiometry[2] = 2.0;  // NH3
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 5331.9;    // K
  reaction_parameters->Keq_coeffs[0] /= SC->TEMP; // non-dimensional
  reaction_parameters->Keq_coeffs[1] = -11.884;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0;

  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    reaction_parameters->volatiles[i] = -1; /* error value */
  for (v = 0; v < Ap->n_volatiles; ++v)
  {
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "N2", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "H2", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[1] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "NH3", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[2] = v;
  }
  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    if (reaction_parameters->volatiles[i] == -1)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Didn't find required volatiles for reaction %s", reaction_parameters->type);
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersCreateWaterIVTANTHERMO(ReactionParameters *reaction_parameters_ptr, const AtmosphereParameters Ap, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  PetscInt i, v;
  PetscBool flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, reaction_parameters_ptr);
  CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "waterIVTANTHERMO";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2, &reaction_parameters->Keq_coeffs, reaction_parameters->n_volatiles, &reaction_parameters->stoichiometry, reaction_parameters->n_volatiles, &reaction_parameters->volatiles);
  CHKERRQ(ierr);

  /* H2O = H2 + 0.5 fO2 */
  /* Keq from IVTANTHERMO for 298.15 < T < 2000 K */
  /* see Schaefer and Fegley (2017) */
  reaction_parameters->stoichiometry[0] = -1.0; // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = -12794;    // K
  reaction_parameters->Keq_coeffs[0] /= SC->TEMP; // non-dimensional
  reaction_parameters->Keq_coeffs[1] = 2.7768;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.5;

  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    reaction_parameters->volatiles[i] = -1; /* error value */
  for (v = 0; v < Ap->n_volatiles; ++v)
  {
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "H2O", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "H2", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[1] = v;
  }
  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    if (reaction_parameters->volatiles[i] == -1)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Didn't find required volatiles for reaction %s", reaction_parameters->type);
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersCreateWaterJANAF(ReactionParameters *reaction_parameters_ptr, const AtmosphereParameters Ap, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  PetscInt i, v;
  PetscBool flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, reaction_parameters_ptr);
  CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "waterJANAF";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2, &reaction_parameters->Keq_coeffs, reaction_parameters->n_volatiles, &reaction_parameters->stoichiometry, reaction_parameters->n_volatiles, &reaction_parameters->volatiles);
  CHKERRQ(ierr);

  /* H2 + 0.5 fO2 = H2O */
  /* Keq from JANAF for 1500 < T < 3000 K */
  reaction_parameters->stoichiometry[0] = -1.0; // H2
  reaction_parameters->stoichiometry[1] = 1.0;  // H2O
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 13152.477779978302; // K, computed from 251801/( (np.log(10)*R))
  reaction_parameters->Keq_coeffs[0] /= SC->TEMP;          // non-dimensional
  reaction_parameters->Keq_coeffs[1] = -3.038586383273608; // computed from -58.173/(np.log(10)*R)
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = -0.5;

  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    reaction_parameters->volatiles[i] = -1; /* error value */
  for (v = 0; v < Ap->n_volatiles; ++v)
  {
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "H2", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "H2O", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[1] = v;
  }
  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    if (reaction_parameters->volatiles[i] == -1)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Didn't find required volatiles for reaction %s", reaction_parameters->type);
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersCreateCarbonDioxideJANAF(ReactionParameters *reaction_parameters_ptr, const AtmosphereParameters Ap, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  PetscInt i, v;
  PetscBool flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, reaction_parameters_ptr);
  CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "carbondioxideJANAF";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2, &reaction_parameters->Keq_coeffs, reaction_parameters->n_volatiles, &reaction_parameters->stoichiometry, reaction_parameters->n_volatiles, &reaction_parameters->volatiles);
  CHKERRQ(ierr);

  /* CO + 0.5 fO2 = CO2 */
  /* Keq from JANAF for 1500 < T < 3000 K */
  reaction_parameters->stoichiometry[0] = -1.0; // CO
  reaction_parameters->stoichiometry[1] = 1.0;  // CO2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 14467.511400133637; // K, computed from 276977/(np.log(10)*R)
  reaction_parameters->Keq_coeffs[0] /= SC->TEMP;          // non-dimensional
  reaction_parameters->Keq_coeffs[1] = -4.348135473316284; // computed from -83.244/(np.log(10)*R)
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = -0.5;

  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    reaction_parameters->volatiles[i] = -1; /* error value */
  for (v = 0; v < Ap->n_volatiles; ++v)
  {
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "CO", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "CO2", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[1] = v;
  }
  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    if (reaction_parameters->volatiles[i] == -1)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Didn't find required volatiles for reaction %s", reaction_parameters->type);
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersCreateCarbonDioxideIVTANTHERMO(ReactionParameters *reaction_parameters_ptr, const AtmosphereParameters Ap, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  PetscInt i, v;
  PetscBool flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, reaction_parameters_ptr);
  CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "carbondioxideIVTANTHERMO";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2, &reaction_parameters->Keq_coeffs, reaction_parameters->n_volatiles, &reaction_parameters->stoichiometry, reaction_parameters->n_volatiles, &reaction_parameters->volatiles);
  CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0; // CO2
  reaction_parameters->stoichiometry[1] = 1.0;  // CO
  /* equilibrium constant coefficients */
  /* Keq from IVTANTHERMO for 298.15 to 2000 K */
  /* see Schaefer and Fegley (2017) */
  reaction_parameters->Keq_coeffs[0] = -14787;    // K
  reaction_parameters->Keq_coeffs[0] /= SC->TEMP; // non-dimensional
  reaction_parameters->Keq_coeffs[1] = 4.5472;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.5;

  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    reaction_parameters->volatiles[i] = -1; /* error value */
  for (v = 0; v < Ap->n_volatiles; ++v)
  {
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "CO2", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v]->prefix, "CO", &flg);
    CHKERRQ(ierr);
    if (flg)
      reaction_parameters->volatiles[1] = v;
  }
  for (i = 0; i < reaction_parameters->n_volatiles; ++i)
    if (reaction_parameters->volatiles[i] == -1)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Didn't find required volatiles for reaction %s", reaction_parameters->type);
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersDestroy(ReactionParameters *reaction_parameters_ptr)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscFree3((*reaction_parameters_ptr)->Keq_coeffs, (*reaction_parameters_ptr)->stoichiometry, (*reaction_parameters_ptr)->volatiles);
  CHKERRQ(ierr);
  ierr = PetscFree(*reaction_parameters_ptr);
  CHKERRQ(ierr); /* must be last */
  *reaction_parameters_ptr = NULL;
  PetscFunctionReturn(0);
}

/* Compute equilibrium constant */
static PetscScalar get_log10_equilibrium_constant(const ReactionParameters reaction_parameters, PetscScalar temp)
{
  /* log10Keq = a/T + b */
  /* e.g., Schaefer and Fegley (2017) */

  PetscScalar log10Keq;

  log10Keq = reaction_parameters->Keq_coeffs[0] / temp + reaction_parameters->Keq_coeffs[1];

  return log10Keq;
}

PetscScalar get_dlog10KdT(const ReactionParameters reaction_parameters, PetscScalar temp)
{
  /* temperature derivation of log10Keq = -a/T^2 */

  PetscScalar dlog10KdT;

  dlog10KdT = -reaction_parameters->Keq_coeffs[0] / PetscPowScalar(temp, 2.0);

  return dlog10KdT;
}

/* Compute modified equilibrium constant */
/* This includes fO2, which helps numerically since the total quantity is better scaled */
PetscScalar get_log10_modified_equilibrium_constant(const ReactionParameters reaction_parameters, PetscScalar temp, const Atmosphere *A, const ScalingConstants SC)
{
  PetscScalar log10G, log10K;

  log10K = get_log10_equilibrium_constant(reaction_parameters, temp);

  log10G = log10K - reaction_parameters->fO2_stoichiometry * A->log10fO2;

  /* account for scaling, since equilibrium constants are defined
     for quantities normalised to 1 bar */
  PetscScalar expon = get_psurf_exponent(reaction_parameters);
  log10G += expon * (5 - PetscLog10Real(SC->PRESSURE));

  return log10G;
}

/* Derivative of log10 modified equilibrium constant with respect to temperature */
PetscScalar get_dlog10GdT(const ReactionParameters reaction_parameters, PetscScalar temp, const Atmosphere *A)
{
  PetscScalar dlog10KdT, dlog10GdT;

  dlog10KdT = get_dlog10KdT(reaction_parameters, A->tsurf);

  dlog10GdT = dlog10KdT - reaction_parameters->fO2_stoichiometry * A->dlog10fO2dT;

  return dlog10GdT;
}

/* Exponent of extra factor of psurf to ensure reaction quotient is scaled correctly */
static PetscScalar get_psurf_exponent(const ReactionParameters reaction_parameters)
{
  PetscInt j;
  PetscScalar expon = 0;

  for (j = 0; j < reaction_parameters->n_volatiles; ++j)
  {
    expon += reaction_parameters->stoichiometry[j];
  }

  return expon;
}

/* Compute reaction quotient (products, numerator) */

PetscScalar get_reaction_quotient_products(const ReactionParameters reaction_parameters, const Atmosphere *A)
{

  return get_reaction_quotient(reaction_parameters, A, 1);
}

PetscScalar get_reaction_quotient_reactants(const ReactionParameters reaction_parameters, const Atmosphere *A)
{

  return get_reaction_quotient(reaction_parameters, A, -1);
}

static PetscScalar get_reaction_quotient(const ReactionParameters reaction_parameters, const Atmosphere *A, PetscInt SIGN)
{
  /* returns numerator for SIGN=1 (products) and denominator for SIGN=-1 (reactants) */
  /* note: excludes fO2 since this is badly scaled.  fO2 is dealt with at the same time
     as the equilibrium constant */

  PetscInt j;
  PetscScalar Q = 1;

  for (j = 0; j < reaction_parameters->n_volatiles; ++j)
  {

    const PetscInt v = reaction_parameters->volatiles[j];

    if (SIGN * reaction_parameters->stoichiometry[j] > 0.0)
    {
      Q *= PetscPowScalar(A->volatiles[v].p, SIGN * reaction_parameters->stoichiometry[j]);
    }
  }

  return Q;
}

PetscScalar get_reaction_quotient_products_time_derivative(const ReactionParameters reaction_parameters, const Atmosphere *A, const AtmosphereParameters Ap)
{

  return get_reaction_quotient_time_derivative(reaction_parameters, A, Ap, 1);
}

PetscScalar get_reaction_quotient_reactants_time_derivative(const ReactionParameters reaction_parameters, const Atmosphere *A, const AtmosphereParameters Ap)
{

  return get_reaction_quotient_time_derivative(reaction_parameters, A, Ap, -1);
}

/* compute reaction quotient (products, numerator) derivative with respect to time t */
static PetscScalar get_reaction_quotient_time_derivative(const ReactionParameters reaction_parameters, const Atmosphere *A, const AtmosphereParameters Ap, PetscInt SIGN)
{
  /* returns dQp/dt for SIGN=1 (products) and dQr/dt for SIGN=-1 (reactants) */

  PetscInt j, k;
  PetscScalar dQdt = 0, dvdt;
  PetscScalar prefactor;

  /* what follows is the chain rule */

  for (j = 0; j < reaction_parameters->n_volatiles; ++j)
  {

    if (SIGN * reaction_parameters->stoichiometry[j] > 0.0)
    {

      /* compute derivative of this volatile */
      const PetscInt v = reaction_parameters->volatiles[j];
      prefactor = 1.0;

      /* find all other volatiles with positive stoichiometry that define the prefactor */
      for (k = 0; k < reaction_parameters->n_volatiles; ++k)
      {
        if (SIGN * reaction_parameters->stoichiometry[k] > 0.0 && k != j)
        {
          const PetscInt v2 = reaction_parameters->volatiles[k];
          prefactor *= PetscPowScalar(A->volatiles[v2].p, SIGN * reaction_parameters->stoichiometry[k]);
        }
      }

      /* compute the derivative of this volatile */
      dvdt = A->volatiles[v].dpdt;

      /* so the contribution of this derivative is as follows */
      dvdt *= SIGN * reaction_parameters->stoichiometry[j] * PetscPowScalar(A->volatiles[v].p, SIGN * reaction_parameters->stoichiometry[j] - 1.0);
      dQdt += prefactor * dvdt;
    }
  }

  return dQdt;
}
