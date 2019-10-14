#include "reaction.h"

/* Note: this could logically be included in parameters.c, but that file was getting crowded */

/* Create a simple reaction, involving 2 volatiles, described by "exchange rate"
   factors gamma0 and gamma1, and an equilibrium ratio factors epsilon0 and epsilon1 */

#if 0
PetscErrorCode ReactionParametersCreateSimple(ReactionParameters* reaction_parameters_ptr,PetscInt v0,PetscInt v1,PetscReal gamma0,PetscReal gamma1,PetscReal epsilon0, PetscReal epsilon1)
{
  PetscErrorCode     ierr;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "simple";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->gamma,2,&reaction_parameters->volatiles,2,&reaction_parameters->epsilon);CHKERRQ(ierr);
  reaction_parameters->gamma[0] = gamma0;
  reaction_parameters->gamma[1] = gamma1;
  reaction_parameters->epsilon[0] = epsilon0;
  reaction_parameters->epsilon[1] = epsilon1;
  reaction_parameters->volatiles[0] = v0;
  reaction_parameters->volatiles[1] = v1;
  PetscFunctionReturn(0);
}
#endif

/* A named reaction */
PetscErrorCode ReactionParametersCreateMethane1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "methane1";
  reaction_parameters->n_volatiles = 4;
  ierr = PetscMalloc2(3,&reaction_parameters->stoichiometry,3,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = 1.0;  // CO2
  reaction_parameters->stoichiometry[1] = 2.0;  // H2
  reaction_parameters->stoichiometry[2] = -1.0; // CH4
  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"CO2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"CH4",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[2] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

/* A named reaction */
PetscErrorCode ReactionParametersCreateWater1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "water1";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,3,&reaction_parameters->stoichiometry,3,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = -12794;
  reaction_parameters->Keq_coeffs[1] = 2.7768;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.5;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2O",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersDestroy(ReactionParameters* reaction_parameters_ptr)
{
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = PetscFree3((*reaction_parameters_ptr)->Keq_coeffs,(*reaction_parameters_ptr)->stoichiometry,(*reaction_parameters_ptr)->volatiles);CHKERRQ(ierr);
  ierr = PetscFree(*reaction_parameters_ptr);CHKERRQ(ierr);/* must be last */
  *reaction_parameters_ptr = NULL;
  PetscFunctionReturn(0);
}

/* Compute equilibrium constant */
PetscScalar get_equilibrium_constant( const ReactionParameters * reaction_parameters_ptr, PetscScalar temp, const Constants *C )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscScalar        log10Keq, Keq;

    temp *= C->TEMP;

    /* log10Keq = a/T + b is a standard form for computing an equilibrium constant */
    /* e.g., Schaefer and Fegley (2017) */
    log10Keq = reaction_parameters->Keq_coeffs[0] / temp + reaction_parameters->Keq_coeffs[1];

    Keq = PetscPowScalar( 10.0, log10Keq );

    return Keq;

}

/* Compute reaction quotient (products, numerator) */
PetscScalar get_reaction_quotient_products( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscInt           j;
    PetscScalar        Qp = 1;

    for (j=0; j<reaction_parameters->n_volatiles; ++j) {

        const PetscInt v = reaction_parameters->volatiles[j];

        if( reaction_parameters->stoichiometry[j] > 0.0 ){
            /* Note: reaction quotient numerator (products), excluding fO2 */
            /* NOTE: introduced scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
            /* TODO: if this works, swap out A->volatiles[v0].p/A->psurf for the volume mixing ratio? */
            Qp *= PetscPowScalar( (A->volatiles[v].p/A->psurf), reaction_parameters->stoichiometry[j] );
        }

    }

    if( reaction_parameters->fO2_stoichiometry > 0.0 ){
        if( A->oxygen_fugacity != 0 ){
            Qp *= PetscPowScalar( A->oxygen_fugacity, reaction_parameters->fO2_stoichiometry );
        }
        else{
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Reaction demands fO2, but OXYGEN_FUGACITY is switched OFF!");
        }
    }

    return Qp;

}

/* Compute reaction quotient (reactants, denominator) */
PetscScalar get_reaction_quotient_reactants( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscInt           j;
    PetscScalar        Qr = 1;

    for (j=0; j<reaction_parameters->n_volatiles; ++j) {

        const PetscInt v = reaction_parameters->volatiles[j];

        if( reaction_parameters->stoichiometry[j] < 0.0 ){
            /* Note: reaction quotient denominator (reactants), excluding fO2 */
            /* NOTE: introduced scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
            /* TODO: if this works, swap out A->volatiles[v0].p/A->psurf for the volume mixing ratio? */
            Qr *= PetscPowScalar( (A->volatiles[v].p/A->psurf), -reaction_parameters->stoichiometry[j] );
        }

    }

    if( reaction_parameters->fO2_stoichiometry < 0.0 ){
        if( A->oxygen_fugacity != 0 ){
            /* Note: negative stoichiometry below, since return denominator and not 1/denominator */
            /* NOTE: introduced scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
            /* TODO: if this works, swap out A->volatiles[v0].p/A->psurf for the volume mixing ratio? */
            Qr *= PetscPowScalar( A->oxygen_fugacity, -reaction_parameters->fO2_stoichiometry );
        }
        else{
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Reaction demands fO2, but OXYGEN_FUGACITY is switched OFF!");
        }
    }

    return Qr;

}
