#include "reaction.h"

static PetscScalar get_psurf_exponent( const ReactionParameters * );
static PetscScalar get_reaction_quotient( const ReactionParameters *, const Atmosphere *, PetscInt );
static PetscScalar get_reaction_quotient_time_derivative( const ReactionParameters *, const Atmosphere *, const AtmosphereParameters *, PetscInt );

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
  reaction_parameters->n_volatiles = 3;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // CO2
  reaction_parameters->stoichiometry[1] = -2.0;  // H2
  reaction_parameters->stoichiometry[2] = 1.0; // CH4
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = -16276;
  reaction_parameters->Keq_coeffs[1] = -5.4738;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 1.0;

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
PetscErrorCode ReactionParametersCreateAmmonia1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "ammonia1";
  reaction_parameters->n_volatiles = 3;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // N2
  reaction_parameters->stoichiometry[1] = -3.0;  // H2
  reaction_parameters->stoichiometry[2] = 2.0; // NH3
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 5331.9;
  reaction_parameters->Keq_coeffs[1] = -11.884;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"N2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"NH3",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[2] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

/* A named reaction */
/* This one is used for testing, since it assumes constant fO2 and K */
/* K = pH2/pH2O = 10**2.0 = 100 */
/* simplified version used in thesis work of Taylor Kutra (2019) */
PetscErrorCode ReactionParametersCreateSimpleWater1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "simplewater1";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 0;
  reaction_parameters->Keq_coeffs[1] = 2.0;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.0;

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

/* A named reaction */
/* This one is used for testing, since it assumes constant fO2 and K */
/* K = pH2/pH2O = 10**1.0 = 10 */
PetscErrorCode ReactionParametersCreateSimpleWater2(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "simplewater2";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 0;
  reaction_parameters->Keq_coeffs[1] = 1.0;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.0;

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

/* A named reaction */
/* This one is used for testing, since it assumes constant fO2 and K */
/* K = pH2/pH2O = 10**0.0 = 1 */
PetscErrorCode ReactionParametersCreateSimpleWater3(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "simplewater3";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 0;
  reaction_parameters->Keq_coeffs[1] = 0.0;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.0;

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
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
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

/* A named reaction */
PetscErrorCode ReactionParametersCreateCarbonDioxide1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "carbondioxide1";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // CO2
  reaction_parameters->stoichiometry[1] = 1.0;  // CO
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = -14787;
  reaction_parameters->Keq_coeffs[1] = 4.5472;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.5;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"CO2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"CO",&flg);CHKERRQ(ierr);
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
  ierr = PetscFree(*reaction_parameters_ptr);CHKERRQ(ierr); /* must be last */
  *reaction_parameters_ptr = NULL;
  PetscFunctionReturn(0);
}

/* TODO: it will probably be helpful for debugging and output to store
   the equilibrium constant and its derivative, but currently we do not
   have an obvious place to store them since they are not parameters
   but rather time-dependent quantities */

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

/* Derivative of equilibrium constant with respect to temperature */
PetscScalar get_equilibrium_constant_temperature_derivative( const ReactionParameters * reaction_parameters_ptr, PetscScalar temp, const Constants *C )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscScalar        dKeqdT;

    temp *= C->TEMP;

    dKeqdT = get_equilibrium_constant( reaction_parameters_ptr, temp, C );
    dKeqdT *= PetscLogReal(10.0);
    dKeqdT *= -reaction_parameters->Keq_coeffs[0] / PetscPowScalar( temp, 2.0 );

    return dKeqdT;

}

/* Exponent of extra factor of psurf to ensure reaction quotient is non-dimensional */
static PetscScalar get_psurf_exponent( const ReactionParameters * reaction_parameters_ptr )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscInt    j;
    PetscScalar expon = 0;

    for( j=0; j<reaction_parameters->n_volatiles; ++j) {
        expon += reaction_parameters->stoichiometry[j];
    }

    return expon;

}

/* Compute reaction quotient (products, numerator) */

PetscScalar get_reaction_quotient_products( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A )
{

    return get_reaction_quotient( reaction_parameters_ptr, A, 1 );

}

PetscScalar get_reaction_quotient_reactants( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A )
{

    return get_reaction_quotient( reaction_parameters_ptr, A, -1 );

}

static PetscScalar get_reaction_quotient( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A, PetscInt SIGN )
{
    /* returns numerator for SIGN=1 (products) and denominator for SIGN=-1 (reactants) */

    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscInt           j;
    PetscScalar        Q = 1;
    PetscScalar        expon;

    for (j=0; j<reaction_parameters->n_volatiles; ++j) {

        const PetscInt v = reaction_parameters->volatiles[j];

        if( SIGN * reaction_parameters->stoichiometry[j] > 0.0 ){
            /* NOTE: introduced scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
            /* TODO: if this works, swap out A->volatiles[v0].p/A->psurf for the volume mixing ratio? */
            Q *= PetscPowScalar( (A->volatiles[v].p/A->psurf), SIGN * reaction_parameters->stoichiometry[j] );
        }

    }

    if( SIGN * reaction_parameters->fO2_stoichiometry > 0.0 ){
        if( A->fO2 != 0 ){
            Q *= PetscPowScalar( A->fO2, SIGN * reaction_parameters->fO2_stoichiometry );
        }
        else{
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Reaction demands fO2, but OXYGEN_FUGACITY is switched OFF!");
        }
    }

    /* since we are using the volume mixing ratio (i.e. scaling by A->psurf),
       we must account for the extra factors of A->psurf to ensure that
       Q=Qp/Qr is non-dimensional */
    expon = get_psurf_exponent( reaction_parameters_ptr );
    if( SIGN * expon > 0.0 ){
        Q *= PetscPowScalar( A->psurf, SIGN * expon );
    }

    return Q;

}

PetscScalar get_reaction_quotient_products_time_derivative( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A, const AtmosphereParameters *Ap )
{

    return get_reaction_quotient_time_derivative( reaction_parameters_ptr, A, Ap, 1 );

}

PetscScalar get_reaction_quotient_reactants_time_derivative( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A, const AtmosphereParameters *Ap )
{

    return  get_reaction_quotient_time_derivative( reaction_parameters_ptr, A, Ap, -1 );

}

/* compute reaction quotient (products, numerator) derivative with respect to time t */
static PetscScalar get_reaction_quotient_time_derivative( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A, const AtmosphereParameters *Ap, PetscInt SIGN )
{
    /* returns dQp/dt for SIGN=1 (products) and dQr/dt for SIGN=-1 (reactants) */

    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscInt           i,j,k;
    PetscScalar        dQdt=0, dpsurfdt, dvdt;
    PetscScalar        expon, prefactor;
    PetscBool          INCLUDE_PSURF = PETSC_FALSE, INCLUDE_FO2 = PETSC_FALSE;

    /* this is a bit ugly, because I decide whether to include the scaling of surface pressure in the
       numerator if the stoichiometry is positive (otherwise it is included in the denominator */
    expon = get_psurf_exponent( reaction_parameters_ptr );
    if( SIGN * expon > 0.0 ){
        INCLUDE_PSURF = PETSC_TRUE;
    }

    /* and similarly for fO2 */
    if( SIGN * reaction_parameters->fO2_stoichiometry > 0.0 ){
        INCLUDE_FO2 = PETSC_TRUE;
    }

    /* TODO: dPsurf/dt is required in various locations, and should probably be computed
       once and stored.  For example, this same calculation is in get_dxdt() */
    dpsurfdt = 0.0;
    /* NOTE: this is over all volatiles, and not just those volatiles in a
       particular reaction */
    for (i=0; i<Ap->n_volatiles; ++i) {
        dpsurfdt += A->volatiles[i].dpdx * A->volatiles[i].dxdt;
    }

    /* what follows is the chain rule */

    /* first, account for volatiles in this reaction that are explicitly tracked
       by the code */
    for (j=0; j<reaction_parameters->n_volatiles; ++j) {

        if( SIGN * reaction_parameters->stoichiometry[j] > 0.0 ){

            /* compute derivative of this volatile */
            const PetscInt v = reaction_parameters->volatiles[j];
            prefactor = 1.0;

            /* find all other volatiles with positive stoichiometry that define the prefactor */
            for (k=0; k<reaction_parameters->n_volatiles; ++k) {
                if ( SIGN * reaction_parameters->stoichiometry[k] > 0.0 && k!=j ){
                    const PetscInt v2 = reaction_parameters->volatiles[k];
                    prefactor *= PetscPowScalar( (A->volatiles[v2].p/A->psurf), SIGN * reaction_parameters->stoichiometry[k] );
                }
            }
            /* should the prefactor include fO2? */
            if( INCLUDE_FO2 ){
                prefactor *= PetscPowScalar( A->fO2, SIGN * reaction_parameters->fO2_stoichiometry );
            }
            /* should the prefactor include Psurf? */
            if( INCLUDE_PSURF ){
                prefactor *= PetscPowScalar( A->psurf, SIGN * expon );
            }
            /* compute the derivative of this volatile */
            /* TODO: this is also a generic formula for the time derivative of the
               volume mixing ratio and appears elsewhere in the code */
            dvdt = (1.0/A->psurf) * A->volatiles[v].dxdt * A->volatiles[v].dpdx;
            dvdt -= (A->volatiles[v].p / PetscPowScalar( A->psurf, 2.0 )) * dpsurfdt;

            /* so the contribution of this derivative is as follows */
            dvdt *= SIGN * reaction_parameters->stoichiometry[j] * PetscPowScalar( (A->volatiles[v].p/A->psurf), SIGN * reaction_parameters->stoichiometry[j] - 1.0 );
            dQdt += prefactor * dvdt;
        }
    }

    /* now account for dfO2/dt */
    if( INCLUDE_FO2 ){

        /* find all other volatiles with positive stoichiometry that define the prefactor */
        prefactor = 1.0;
        for (k=0; k<reaction_parameters->n_volatiles; ++k) {
            if ( SIGN * reaction_parameters->stoichiometry[k] > 0.0 ){
                const PetscInt v2 = reaction_parameters->volatiles[k];
                prefactor *= PetscPowScalar( (A->volatiles[v2].p/A->psurf), SIGN * reaction_parameters->stoichiometry[k] );
            }
        }

        /* should the prefactor include Psurf? */
        if( INCLUDE_PSURF ){
            prefactor *= PetscPowScalar( A->psurf, SIGN * expon );
        }

        /* now compute the time derivative of fO2 */
        /* TODO: check that these entries are set and updated accordingly */
        dvdt = A->dfO2dT * A->dtsurfdt;

        /* contribution is */
        dvdt *= SIGN * reaction_parameters->fO2_stoichiometry * PetscPowScalar( A->fO2, SIGN * reaction_parameters->fO2_stoichiometry - 1.0 );
        dQdt += prefactor * dvdt;
    }

    /* now account for dpsurf/dt */
    if( INCLUDE_PSURF ) {

        /* find all other volatiles with positive stoichiometry that define the prefactor */
        prefactor = 1.0;
        for (k=0; k<reaction_parameters->n_volatiles; ++k) {
            if ( SIGN * reaction_parameters->stoichiometry[k] > 0.0 ){
                const PetscInt v2 = reaction_parameters->volatiles[k];
                prefactor *= PetscPowScalar( (A->volatiles[v2].p/A->psurf), SIGN * reaction_parameters->stoichiometry[k] );
            }
        }

        /* should the prefactor include fO2? */
        if( INCLUDE_FO2 ){
            prefactor *= PetscPowScalar( A->fO2, SIGN * reaction_parameters->fO2_stoichiometry );
        }

        /* now compute the time derivative of psurf */
        dvdt = dpsurfdt;

        /* so the contribution is */
        dvdt *= SIGN * expon * PetscPowScalar( A->psurf, SIGN * expon - 1.0 );
        dQdt += prefactor * dvdt;
    }

    return dQdt;

}
